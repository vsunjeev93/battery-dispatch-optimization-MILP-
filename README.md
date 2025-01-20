# Battery Energy Dispatch Optimizer

Optimizes battery storage and grid interaction for homes with solar panels to minimize electricity costs while participating in demand response programs.

## Requirements

- Python 3.x
- Dependencies:
  - matplotlib
  - mip
  - pandas
  - datetime

## Input Files

### input_profiles.csv
CSV file with 48 rows (half-hourly data) containing:
- Time: Timestamp (format: "MM/DD/YY HH:MM")
- Load (kW): Household power demand
- PV (kW): Solar panel generation

## Configuration Parameters

### Battery Parameters
- Charge/Discharge Efficiency: 95%
- Maximum Charge/Discharge Power: ±25 kW
- Maximum Energy Capacity: 53 kWh
- Initial Energy: 0 kWh

### Time Windows
- Peak Demand: 17:00-21:00 (5 PM-9 PM)
- Demand Response: 19:00-20:00 (7 PM-8 PM)

### Cost Parameters
- Import Cost: $0.10/kWh
- Export Revenue: $0.03/kWh
- Peak Demand Charge: $9/kW
- Demand Response Revenue: $10/kWh

## Usage

```python
python dispatch.py
```

## Outputs

1. results.csv: Optimized operation schedule
2. results.png: Visualization showing:
   - Grid power flow
   - Battery power
   - Battery energy level

## Optimization Features

- Minimizes electricity costs
- Manages battery charging/discharging
- Optimizes solar power usage
- Participates in demand response events
- Avoids peak demand charges
- Prevents simultaneous charging/discharging
