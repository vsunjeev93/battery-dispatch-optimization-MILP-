#!/usr/bin/python
import matplotlib.pyplot as plt
from mip import Model, xsum, BINARY, CONTINUOUS, minimize, OptimizationStatus
from datetime import datetime
import pandas as pd


class OptimizeDispatch:
   def __init__(
       self,
       profiles,
       charge_eff,
       discharge_eff,
       max_charge,
       max_discharge,
       max_energy,
       initial_energy,
       peak_demand_begin,
       peak_demand_end,
       demand_response_begin,
       demand_response_end,
   ):
       # params
       self.profiles = profiles
       self.charge_eff = charge_eff
       self.discharge_eff = discharge_eff
       self.max_charge = max_charge
       self.max_discharge = max_discharge
       self.max_energy = max_energy
       self.initial_energy = initial_energy
       self.peak_demand_begin = peak_demand_begin
       self.peak_demand_end = peak_demand_end
       self.demand_response_begin = demand_response_begin
       self.demand_response_end = demand_response_end

   def get_time_difference(self, time1, time2):
       """helper method to get difference between two time strings in hrs"""
       datetime_object1 = datetime.strptime(time1, "%m/%d/%y %H:%M")
       datetime_object2 = datetime.strptime(time2, "%m/%d/%y %H:%M")
       difference = datetime_object2 - datetime_object1
       return float(difference.total_seconds() / 3600)

   def parse_profiles(self):
       """method to read the input csv"""
       profiles_df = pd.read_csv(self.profiles)
       household_load = {}
       pv_gen_power = {}
       time_stamps = []
       for index, row in profiles_df.iterrows():
           household_load[row["Time"]] = row["Load (kW)"]
           pv_gen_power[row["Time"]] = row["PV (kW)"]
           time_stamps.append(row["Time"])

       return household_load, pv_gen_power, time_stamps

   def run_optimization(self, c_import, r_export, c_peak_demand, r_demand_response):
       """defines decisions variables, constraints and runs the optimization"""
       h_load, P_pv, times = self.parse_profiles()
       solver = Model()
       M1 = max(h_load.values()) + max(P_pv.values())
       M2 = max(self.max_charge, -self.max_discharge)
       dt = self.get_time_difference(
           times[0], times[1]
       )  # getting time between two loads (assumption is that this is constant)

       # decisison variables
       Pb = {}  # power from battery -ve is discharging +ve is charging
       P_grid = {}  # power from grid -ve is export +ve is import
       Eb = {}  # energy stored in battery
       Pb_in = {}  # auxillary variable for battery charging always >=0
       Pb_out = {}  # auxillary variable for  battery discharging >=0
       P_import = {}  # auxillary variable for power imported from grid >=0
       P_export = {}  # auxillary variable for power exported from grid >=0

       """decision variables for big M constraints
       (P_import/P_export cannot be simultaneously nonzero)"""
       delta = {}
       beta = {}

       """ decision variables for big M constraints
       (Pb_in/Pb_out cannot be simultaneously nonzero)"""
       theta = {}
       gamma = {}

       for t, P_load in h_load.items():
           # creating the variables
           Pb[t] = solver.add_var(
               f"Pb_{t}",
               var_type=CONTINUOUS,
               lb=self.max_discharge,
               ub=self.max_charge,
           )
           Eb[t] = solver.add_var(f"Eb_{t}", var_type=CONTINUOUS, ub=self.max_energy)
           Pb_in[t] = solver.add_var(
               f"Pb_in_{t}", var_type=CONTINUOUS, ub=self.max_charge
           )
           Pb_out[t] = solver.add_var(
               f"Pb_out_{t}", var_type=CONTINUOUS, ub=-self.max_discharge
           )
           P_grid[t] = solver.add_var(
               f"P_grid_{t}", var_type=CONTINUOUS, lb=-float("inf")
           )
           P_import[t] = solver.add_var(f"P_import_{t}", var_type=CONTINUOUS)
           P_export[t] = solver.add_var(f"P_export_{t}", var_type=CONTINUOUS)
           # P_grid=P_import-P_export
           solver += P_import[t] - P_export[t] == P_grid[t]
           # P_grid= sum of load plus battery and pv power
           solver += P_load + Pb[t] - P_pv[t] == P_grid[t]
           # Pb_in-Pb_out==Pb (net)
           solver += Pb_in[t] - Pb_out[t] == Pb[t]
           # cannot charge more than avaiable PV power
           solver += Pb[t] <= P_pv[t]
       # a variable to capture the max power in the peak demand phase
       z = solver.add_var("z", var_type=CONTINUOUS)
       # defining a terminal battery state at midnight next day
       Eb_terminal = solver.add_var("Eb_T", var_type=CONTINUOUS, ub=self.max_energy)
       for index, t in enumerate(times):
           if index == 0:
               # inital battery energy=0
               solver += Eb[t] == self.initial_energy
           if index == len(times) - 1:
               # battery dynamics for the last time point
               solver += (
                   Eb_terminal
                   == Eb[t]
                   + Pb_in[t] * (self.charge_eff) * dt
                   - Pb_out[t] / (self.discharge_eff) * dt
               )
           else:
               # dynamics for all other times
               solver += (
                   Eb[times[index + 1]]
                   == Eb[t]
                   + Pb_in[t] * (self.charge_eff) * dt
                   - Pb_out[t] / (self.discharge_eff) * dt
               )
           # getting the max of power imported in the peak demand phase
           if (
               self.get_time_difference(times[0], t) >= self.peak_demand_begin
               and self.get_time_difference(times[0], t) < self.peak_demand_end
           ):
               solver += z >= P_import[t]
       """'simultaneous import/export from/to grid not possible
       ( problem is unbounded if we don't consider this)"""
       for t in times:
           delta[t] = solver.add_var(f"delta_{t}", var_type=BINARY)
           beta[t] = solver.add_var(f"beta_{t}", var_type=BINARY)
           solver += P_import[t] <= M1 * delta[t]
           solver += P_export[t] <= M1 * beta[t]
           solver += beta[t] + delta[t] <= 1
       # simultaneous charing/discharging is not physically possible
       for t in times:
           theta[t] = solver.add_var(f"theta_{t}", var_type=BINARY)
           gamma[t] = solver.add_var(f"gamma_{t}", var_type=BINARY)
           solver += Pb_in[t] <= M2 * theta[t]
           solver += Pb_out[t] <= M2 * gamma[t]
           solver += gamma[t] + theta[t] <= 1
       # peak demand cost contribution to objective
       obj = [c_peak_demand * z]
       # import cost
       obj += [c_import * P_import[t] * dt for t in times]
       # export revenue
       obj += [-r_export * P_export[t] * dt for t in times]
       # demand response revenue
       obj += [
           -r_demand_response * P_export[t] * dt
           for t in times
           if self.get_time_difference(times[0], t) >= self.demand_response_begin
           and self.get_time_difference(times[0], t) < self.demand_response_end
       ]
       solver.objective = minimize(xsum(obj))
       status = solver.optimize()
       if (
           status == OptimizationStatus.OPTIMAL
           or status == OptimizationStatus.FEASIBLE
       ):
           print(f"MIP finished {status} {solver.objective_value}")
           self.plot_and_write_results(Pb, P_grid, Eb, times)
       else:
           print(f"No Solution found: {status}")

   def plot_and_write_results(self, battery_power, meter, battery_energy, times):
       """method to plot and write results"""
       battery_power_list = [battery_power[t].x for t in times]
       meter_list = [float(meter[t].x) for t in times]
       battery_energy_list = [float(battery_energy[t].x) for t in times]
       time_in_hrs = [self.get_time_difference(times[0], t) for t in times]
       results_dict = {
           "Time": time_in_hrs,
           "battery_power (kW)": battery_power_list,
           "meter (kW)": meter_list,
       }
       df = pd.DataFrame.from_dict(results_dict, orient="index").transpose()
       df.to_csv("results.csv", index=False)
       plt.figure()
       plt.xlabel("time (hr)", fontsize=20)
       plt.plot(time_in_hrs, meter_list, color="black", label="meter (kW)")
       plt.plot(
           time_in_hrs, battery_power_list, color="green", label="battery_power (kW)"
       )
       plt.plot(
           time_in_hrs, battery_energy_list, color="red", label="battery_energy (kWh)"
       )
       plt.xticks(range(0, 25))
       plt.legend(prop={"size": 15})
       plt.savefig("results.png")
       plt.show()


if __name__ == "__main__":
   cost_params = {
       "c_import": 0.1,
       "r_export": 0.03,
       "c_peak_demand": 9,
       "r_demand_response": 10,
   }
   input_params = {
       "profiles": "input_profiles.csv",
       "charge_eff": 0.95,  # charge efficiency
       "discharge_eff": 0.95,  # discharge efficiency
       "max_charge": 25,  # max charge power of battery
       "max_discharge": -25,  # max discharge power of battery
       "max_energy": 53,  # max energy stored in battery
       "initial_energy": 0,  # initial energy of battery
       "peak_demand_begin": 17,  # 24 hr time notation
       "peak_demand_end": 21,
       "demand_response_begin": 19,
       "demand_response_end": 20,
   }
   dispatcher = OptimizeDispatch(**input_params)
   dispatcher.run_optimization(**cost_params)