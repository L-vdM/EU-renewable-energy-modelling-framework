import pandas as pd
import pyomo.environ as pyo

# Create a solver
opt = pyo.SolverFactory("glpk")
import config


def _model_to_dataframe(model):
    df = pd.DataFrame()
    df["Eout"] = model.Eout.extract_values().values()
    df["netto_demand"] = model.Demand.extract_values().values()
    df["Ein"] = model.Ein.extract_values().values()
    df["inflow"] = model.inflow.extract_values().values()
    df["reservoir"] = model.Res.extract_values().values()
    df["conventional"] = model.Conv.extract_values().values()
    df["diff_dev"] = model.diff_dev.extract_values().values()
    return df


def optimize_reservoirs_usage(
    df,
    cap_hydro,
    cap_conv=None,
    stor_weeks=None,
    variables=["netto_demand", "inflow"],
    begin_fraction=0.0,
    ts=24,
    cyclus=None,
    cost_hydro=1,
    cost_conv=10,
):
    """
    Optimize the charge/discharge of reservoirs per country based on
    reservoir inflow and demand.

    Parameters
    ----------
    df : dataframe
        dataframe with columns of 'demand' and 'inflow'
    hyro_cap: hydropower capacity [MW]
    conv_cap: conventional installed capacity [MWh]


    Returns
    -------
    dataframe
        hourly state of charge, charge/discharge behavior, lbmp, and time stamp
    """
    # reset index
    # remove na values
    df = df[variables].dropna()
    meandemand = df[variables[0]].mean()
    meaninflow = df[variables[1]].mean()
    dft = df.reset_index()
    if cyclus == None:
        cyclus = (0, dft.shape[0])
    model = pyo.ConcreteModel()
    ### set parameters
    model.T = pyo.Set(doc="timestep", initialize=dft.index.tolist(), ordered=True)
    model.CostHydro = pyo.Param(initialize=cost_hydro, doc="cost of hydropower")
    model.CostConv = pyo.Param(initialize=cost_conv, doc="cost of fossil fuels")
    model.CapHydro = pyo.Param(
        initialize=cap_hydro * ts, doc="Max rate of hydropower flow (MWh/day) out"
    )
    if cap_conv is not None:
        model.CapConv = pyo.Param(
            initialize=cap_conv * ts, doc="Max rate of hydropower flow (MWh/day) out"
        )
    if stor_weeks is None:
        # estimate reservoir capacity based on mean 3 week energy demand
        ### TODO: check how much in normal here
        cap_reservoir = df[[variables[1]]].resample("1Y").sum().mean()[0] / 12
    else:
        cap_reservoir = df[[variables[1]]].resample("1Y").sum().mean()[0] / (
            52 / stor_weeks
        )

    model.ResMax = pyo.Param(initialize=cap_reservoir, doc="Max storage (MWh)")
    model.Demand = pyo.Param(
        model.T, initialize=dft[variables[0]].to_list(), doc="daily energy demand"
    )
    model.inflow = pyo.Param(
        model.T, initialize=dft[variables[1]].to_list(), doc="hydropower daily inflow"
    )

    ### set variables
    # reservoirs Ein, Eout, fossilfuel production and storage
    model.Ein = pyo.Var(model.T, domain=pyo.NonNegativeReals)
    model.Eout = pyo.Var(model.T, bounds=(0, model.CapHydro))
    model.Res = pyo.Var(model.T, bounds=(0, model.ResMax))  # storage MWh
    model.demand_deviation = pyo.Var(model.T, domain=pyo.Reals)
    model.flow_deviation = pyo.Var(model.T, domain=pyo.Reals)
    model.diff_dev = pyo.Var(model.T, domain=pyo.NonNegativeReals)
    model.ConvMax = pyo.Var(domain=pyo.NonNegativeReals)

    if cap_conv is not None:
        model.Conv = pyo.Var(model.T, bounds=(0, model.CapConv))
    else:
        model.Conv = pyo.Var(model.T, domain=pyo.NonNegativeReals)

    ### define constraints
    def charge_constraint(model, t):
        "Inflow into reservoir not larger than hydrologic inflow"
        return model.Ein[t] <= model.inflow[t]

    def storage_state(model, t):
        "Amount of energy stored in reservoir"
        # Set first hour state of charge to defined fraction of max
        # note that efficiency is already included in inflow calculations
        if t == model.T.first():
            return model.Res[t] == model.ResMax * begin_fraction
        else:
            return model.Res[t] == (
                model.Res[t - 1] + model.Ein[t - 1] - model.Eout[t - 1]
            )

    def final_storage_state(model, t):
        "Amount of energy stored in reservoir at the end of timeseries"
        # set as the same as the begnning of timeseries
        if t == model.T.last():
            return model.Res[t] >= (model.ResMax * begin_fraction)
        else:
            return pyo.Constraint.Skip

    def positive_charge(model, t):
        "Limit discharge to the amount of charge in reservoir"
        return model.Eout[t] <= model.Res[t]

    def meet_demand(model, t):
        "Meet the demand hydro and conventional energy"
        return model.Demand[t] <= (model.Eout[t] + model.Conv[t])

    #     def yearly_storage_state(model, t):
    #         "Set storage state to be the same at the end of cyclus as beginning"
    #         if (t-cyclus[0]) % cyclus[1] == 0:
    #             return model.Res[t] <= model.ResMax * begin_fraction
    #         else:
    #             return pyo.Constraint.Skip

    def Peak_Rule(model, t):
        return model.ConvMax >= model.Conv[t]

    ## linear expression for absolute value
    def demand_deviation(model, t):
        return model.demand_deviation[t] == (model.Demand[t] - meandemand) / meandemand

    def flow_deviation(model, t):
        return model.flow_deviation[t] == (model.Eout[t] - meaninflow) / meaninflow

    def diff_deviation1(model, t):
        return model.diff_dev[t] >= (
            model.demand_deviation[t] - model.flow_deviation[t]
        )

    def diff_deviation2(model, t):
        return model.diff_dev[t] >= -(
            model.demand_deviation[t] - model.flow_deviation[t]
        )

    ### set constraints
    model.charge = pyo.Constraint(model.T, rule=charge_constraint)
    model.stor_state = pyo.Constraint(model.T, rule=storage_state)
    model.stor_end_state = pyo.Constraint(model.T, rule=final_storage_state)
    model.positive_charge = pyo.Constraint(model.T, rule=positive_charge)
    model.production = pyo.Constraint(model.T, rule=meet_demand)
    #     model.yearly_cycle = pyo.Constraint(model.T, rule=yearly_storage_state)

    model.peak = pyo.Constraint(model.T, rule=Peak_Rule)

    model.demand_dev = pyo.Constraint(model.T, rule=demand_deviation)
    model.flow_dev = pyo.Constraint(model.T, rule=flow_deviation)
    model.diff_dev1 = pyo.Constraint(model.T, rule=diff_deviation1)
    model.diff_dev2 = pyo.Constraint(model.T, rule=diff_deviation2)
    #     model.demand_dev2 = pyo.Constraint(model.T, rule=demand_deviation2)
    #     model.flow_dev1 = pyo.Constraint(model.T, rule=flow_deviation1)
    #     model.flow_dev2 = pyo.Constraint(model.T, rule=flow_deviation2)

    #     ### set cost function
    #     cost = sum((model.CostConv * model.Conv[t] + model.CostHydro * model.Eout[t] + model.diff_dev[t]*meaninflow) for t in model.T)

    cost = (
        sum(
            (
                (model.CostConv * model.Conv[t] + model.CostHydro * model.Eout[t])
                + model.diff_dev[t]
            )
            for t in model.T
        )
    ) + model.ConvMax
    ### set cost function
    #     cost = (sum((model.CostConv * model.Conv[t] + model.CostHydro * model.Eout[t]) for t in model.T))
    #      + sum(model.diff_dev[t] for t in model.T) + model.ConvMax
    ### solve model
    model.objective = pyo.Objective(expr=cost, sense=pyo.minimize)
    opt.solve(model)

    ### save results in dataframe and set index back to timesteps.
    df_results = _model_to_dataframe(model)
    df_results.index = df.index

    return model, df_results
