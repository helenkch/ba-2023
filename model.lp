\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 y[(1,_0)] + 4 y[(3,_0)] + 4 y[(5,_0)] + 8 y[(6,_0)] + 4 y[(1,_1)]
   + 4 y[(3,_1)] + 4 y[(5,_1)] + 8 y[(6,_1)] + 4 y[(1,_2)] + 4 y[(3,_2)]
   + 4 y[(5,_2)] + 8 y[(6,_2)] + 4 y[(3,_3)] + 4 y[(5,_3)] + 8 y[(6,_3)]
   + 4 y[(1,_4)] + 4 y[(2,_4)] + 4 y[(4,_4)] + 8 y[(5,_4)] + 12 y[(6,_4)]
   + 0.8 absDeviationColSum[0] + 0.8 absDeviationColSum[1]
   + 0.8 absDeviationColSum[2] + 0.8 absDeviationColSum[3]
   + 0.8 absDeviationColSum[4] + 0.8 absDeviationColSum[5]
   + 0.8 absDeviationColSum[6] + 0.7 absDeviation[(0,_0)]
   + 0.7 absDeviation[(1,_0)] + 0.7 absDeviation[(2,_0)]
   + 0.7 absDeviation[(3,_0)] + 0.7 absDeviation[(4,_0)]
   + 0.7 absDeviation[(5,_0)] + 0.7 absDeviation[(6,_0)]
   + 0.7 absDeviation[(0,_1)] + 0.7 absDeviation[(1,_1)]
   + 0.7 absDeviation[(2,_1)] + 0.7 absDeviation[(3,_1)]
   + 0.7 absDeviation[(4,_1)] + 0.7 absDeviation[(5,_1)]
   + 0.7 absDeviation[(6,_1)] + 0.7 absDeviation[(0,_2)]
   + 0.7 absDeviation[(1,_2)] + 0.7 absDeviation[(2,_2)]
   + 0.7 absDeviation[(3,_2)] + 0.7 absDeviation[(4,_2)]
   + 0.7 absDeviation[(5,_2)] + 0.7 absDeviation[(6,_2)]
   + 0.7 absDeviation[(0,_3)] + 0.7 absDeviation[(1,_3)]
   + 0.7 absDeviation[(2,_3)] + 0.7 absDeviation[(3,_3)]
   + 0.7 absDeviation[(4,_3)] + 0.7 absDeviation[(5,_3)]
   + 0.7 absDeviation[(6,_3)] + 0.7 absDeviation[(0,_4)]
   + 0.7 absDeviation[(1,_4)] + 0.7 absDeviation[(2,_4)]
   + 0.7 absDeviation[(3,_4)] + 0.7 absDeviation[(4,_4)]
   + 0.7 absDeviation[(5,_4)] + 0.7 absDeviation[(6,_4)] - 71.8474 Constant
Subject To
 min_constraint[(0,_0)]: x[(0,_0)] - 2 y[(0,_0)] >= 0
 max_constraint[(0,_0)]: x[(0,_0)] - 4 y[(0,_0)] <= 0
 min_constraint[(1,_0)]: x[(1,_0)] - 2 y[(1,_0)] >= 0
 max_constraint[(1,_0)]: x[(1,_0)] - 4 y[(1,_0)] <= 0
 min_constraint[(2,_0)]: x[(2,_0)] - 2 y[(2,_0)] >= 0
 max_constraint[(2,_0)]: x[(2,_0)] - 4 y[(2,_0)] <= 0
 min_constraint[(3,_0)]: x[(3,_0)] - 2 y[(3,_0)] >= 0
 max_constraint[(3,_0)]: x[(3,_0)] - 4 y[(3,_0)] <= 0
 min_constraint[(4,_0)]: x[(4,_0)] - 2 y[(4,_0)] >= 0
 max_constraint[(4,_0)]: x[(4,_0)] - 4 y[(4,_0)] <= 0
 min_constraint[(5,_0)]: x[(5,_0)] - 2 y[(5,_0)] >= 0
 max_constraint[(5,_0)]: x[(5,_0)] - 4 y[(5,_0)] <= 0
 min_constraint[(6,_0)]: x[(6,_0)] - 2 y[(6,_0)] >= 0
 max_constraint[(6,_0)]: x[(6,_0)] - 4 y[(6,_0)] <= 0
 min_constraint[(0,_1)]: x[(0,_1)] - 2 y[(0,_1)] >= 0
 max_constraint[(0,_1)]: x[(0,_1)] - 4 y[(0,_1)] <= 0
 min_constraint[(1,_1)]: x[(1,_1)] - 2 y[(1,_1)] >= 0
 max_constraint[(1,_1)]: x[(1,_1)] - 4 y[(1,_1)] <= 0
 min_constraint[(2,_1)]: x[(2,_1)] - 2 y[(2,_1)] >= 0
 max_constraint[(2,_1)]: x[(2,_1)] - 4 y[(2,_1)] <= 0
 min_constraint[(3,_1)]: x[(3,_1)] - 2 y[(3,_1)] >= 0
 max_constraint[(3,_1)]: x[(3,_1)] - 4 y[(3,_1)] <= 0
 min_constraint[(4,_1)]: x[(4,_1)] - 2 y[(4,_1)] >= 0
 max_constraint[(4,_1)]: x[(4,_1)] - 4 y[(4,_1)] <= 0
 min_constraint[(5,_1)]: x[(5,_1)] - 2 y[(5,_1)] >= 0
 max_constraint[(5,_1)]: x[(5,_1)] - 4 y[(5,_1)] <= 0
 min_constraint[(6,_1)]: x[(6,_1)] - 2 y[(6,_1)] >= 0
 max_constraint[(6,_1)]: x[(6,_1)] - 4 y[(6,_1)] <= 0
 min_constraint[(0,_2)]: x[(0,_2)] - 2 y[(0,_2)] >= 0
 max_constraint[(0,_2)]: x[(0,_2)] - 4 y[(0,_2)] <= 0
 min_constraint[(1,_2)]: x[(1,_2)] - 2 y[(1,_2)] >= 0
 max_constraint[(1,_2)]: x[(1,_2)] - 4 y[(1,_2)] <= 0
 min_constraint[(2,_2)]: x[(2,_2)] - 2 y[(2,_2)] >= 0
 max_constraint[(2,_2)]: x[(2,_2)] - 4 y[(2,_2)] <= 0
 min_constraint[(3,_2)]: x[(3,_2)] - 2 y[(3,_2)] >= 0
 max_constraint[(3,_2)]: x[(3,_2)] - 4 y[(3,_2)] <= 0
 min_constraint[(4,_2)]: x[(4,_2)] - 2 y[(4,_2)] >= 0
 max_constraint[(4,_2)]: x[(4,_2)] - 4 y[(4,_2)] <= 0
 min_constraint[(5,_2)]: x[(5,_2)] - 2 y[(5,_2)] >= 0
 max_constraint[(5,_2)]: x[(5,_2)] - 4 y[(5,_2)] <= 0
 min_constraint[(6,_2)]: x[(6,_2)] - 2 y[(6,_2)] >= 0
 max_constraint[(6,_2)]: x[(6,_2)] - 4 y[(6,_2)] <= 0
 min_constraint[(0,_3)]: x[(0,_3)] - 2 y[(0,_3)] >= 0
 max_constraint[(0,_3)]: x[(0,_3)] - 4 y[(0,_3)] <= 0
 min_constraint[(1,_3)]: x[(1,_3)] - 2 y[(1,_3)] >= 0
 max_constraint[(1,_3)]: x[(1,_3)] - 4 y[(1,_3)] <= 0
 min_constraint[(2,_3)]: x[(2,_3)] - 2 y[(2,_3)] >= 0
 max_constraint[(2,_3)]: x[(2,_3)] - 4 y[(2,_3)] <= 0
 min_constraint[(3,_3)]: x[(3,_3)] - 2 y[(3,_3)] >= 0
 max_constraint[(3,_3)]: x[(3,_3)] - 4 y[(3,_3)] <= 0
 min_constraint[(4,_3)]: x[(4,_3)] - 2 y[(4,_3)] >= 0
 max_constraint[(4,_3)]: x[(4,_3)] - 4 y[(4,_3)] <= 0
 min_constraint[(5,_3)]: x[(5,_3)] - 2 y[(5,_3)] >= 0
 max_constraint[(5,_3)]: x[(5,_3)] - 4 y[(5,_3)] <= 0
 min_constraint[(6,_3)]: x[(6,_3)] - 2 y[(6,_3)] >= 0
 max_constraint[(6,_3)]: x[(6,_3)] - 4 y[(6,_3)] <= 0
 min_constraint[(0,_4)]: x[(0,_4)] - 2 y[(0,_4)] >= 0
 max_constraint[(0,_4)]: x[(0,_4)] - 4 y[(0,_4)] <= 0
 min_constraint[(1,_4)]: x[(1,_4)] - 2 y[(1,_4)] >= 0
 max_constraint[(1,_4)]: x[(1,_4)] - 4 y[(1,_4)] <= 0
 min_constraint[(2,_4)]: x[(2,_4)] - 2 y[(2,_4)] >= 0
 max_constraint[(2,_4)]: x[(2,_4)] - 4 y[(2,_4)] <= 0
 min_constraint[(3,_4)]: x[(3,_4)] - 2 y[(3,_4)] >= 0
 max_constraint[(3,_4)]: x[(3,_4)] - 4 y[(3,_4)] <= 0
 min_constraint[(4,_4)]: x[(4,_4)] - 2 y[(4,_4)] >= 0
 max_constraint[(4,_4)]: x[(4,_4)] - 4 y[(4,_4)] <= 0
 min_constraint[(5,_4)]: x[(5,_4)] - 2 y[(5,_4)] >= 0
 max_constraint[(5,_4)]: x[(5,_4)] - 4 y[(5,_4)] <= 0
 min_constraint[(6,_4)]: x[(6,_4)] - 2 y[(6,_4)] >= 0
 max_constraint[(6,_4)]: x[(6,_4)] - 4 y[(6,_4)] <= 0
 capacity_constraint_tours[0]: x[(0,_0)] + x[(1,_0)] + x[(2,_0)]
   + x[(3,_0)] + x[(4,_0)] + x[(5,_0)] + x[(6,_0)] = 12
 capacity_constraint_days[0]: y[(0,_0)] + y[(1,_0)] + y[(2,_0)] + y[(3,_0)]
   + y[(4,_0)] + y[(5,_0)] + y[(6,_0)] = 3
 capacity_constraint_tours[1]: x[(0,_1)] + x[(1,_1)] + x[(2,_1)]
   + x[(3,_1)] + x[(4,_1)] + x[(5,_1)] + x[(6,_1)] = 10
 capacity_constraint_days[1]: y[(0,_1)] + y[(1,_1)] + y[(2,_1)] + y[(3,_1)]
   + y[(4,_1)] + y[(5,_1)] + y[(6,_1)] = 3
 capacity_constraint_tours[2]: x[(0,_2)] + x[(1,_2)] + x[(2,_2)]
   + x[(3,_2)] + x[(4,_2)] + x[(5,_2)] + x[(6,_2)] = 12
 capacity_constraint_days[2]: y[(0,_2)] + y[(1,_2)] + y[(2,_2)] + y[(3,_2)]
   + y[(4,_2)] + y[(5,_2)] + y[(6,_2)] = 3
 capacity_constraint_tours[3]: x[(0,_3)] + x[(1,_3)] + x[(2,_3)]
   + x[(3,_3)] + x[(4,_3)] + x[(5,_3)] + x[(6,_3)] = 15
 capacity_constraint_days[3]: y[(0,_3)] + y[(1,_3)] + y[(2,_3)] + y[(3,_3)]
   + y[(4,_3)] + y[(5,_3)] + y[(6,_3)] = 4
 capacity_constraint_tours[4]: x[(0,_4)] + x[(1,_4)] + x[(2,_4)]
   + x[(3,_4)] + x[(4,_4)] + x[(5,_4)] + x[(6,_4)] = 8
 capacity_constraint_days[4]: y[(0,_4)] + y[(1,_4)] + y[(2,_4)] + y[(3,_4)]
   + y[(4,_4)] + y[(5,_4)] + y[(6,_4)] = 2
 we/publVac_constraint[5]: y[(5,_0)] = 0
 we/publVac_constraint[5]: y[(5,_1)] = 0
 we/publVac_constraint[5]: y[(5,_2)] = 0
 we/publVac_constraint[5]: y[(5,_3)] = 0
 we/publVac_constraint[5]: y[(5,_4)] = 0
 we/publVac_constraint[6]: y[(6,_0)] = 0
 we/publVac_constraint[6]: y[(6,_1)] = 0
 we/publVac_constraint[6]: y[(6,_2)] = 0
 we/publVac_constraint[6]: y[(6,_3)] = 0
 we/publVac_constraint[6]: y[(6,_4)] = 0
 oneVehicle_constraint[0]: y[(0,_0)] <= 1
 oneVehicle_constraint[1]: y[(0,_1)] + y[(0,_4)] <= 1
 oneVehicle_constraint[2]: y[(0,_2)] <= 1
 oneVehicle_constraint[3]: y[(0,_3)] <= 1
 oneVehicle_constraint[0]: y[(1,_0)] <= 1
 oneVehicle_constraint[1]: y[(1,_1)] + y[(1,_4)] <= 1
 oneVehicle_constraint[2]: y[(1,_2)] <= 1
 oneVehicle_constraint[3]: y[(1,_3)] <= 1
 oneVehicle_constraint[0]: y[(2,_0)] <= 1
 oneVehicle_constraint[1]: y[(2,_1)] + y[(2,_4)] <= 1
 oneVehicle_constraint[2]: y[(2,_2)] <= 1
 oneVehicle_constraint[3]: y[(2,_3)] <= 1
 oneVehicle_constraint[0]: y[(3,_0)] <= 1
 oneVehicle_constraint[1]: y[(3,_1)] + y[(3,_4)] <= 1
 oneVehicle_constraint[2]: y[(3,_2)] <= 1
 oneVehicle_constraint[3]: y[(3,_3)] <= 1
 oneVehicle_constraint[0]: y[(4,_0)] <= 1
 oneVehicle_constraint[1]: y[(4,_1)] + y[(4,_4)] <= 1
 oneVehicle_constraint[2]: y[(4,_2)] <= 1
 oneVehicle_constraint[3]: y[(4,_3)] <= 1
 oneVehicle_constraint[0]: y[(5,_0)] <= 1
 oneVehicle_constraint[1]: y[(5,_1)] + y[(5,_4)] <= 1
 oneVehicle_constraint[2]: y[(5,_2)] <= 1
 oneVehicle_constraint[3]: y[(5,_3)] <= 1
 oneVehicle_constraint[0]: y[(6,_0)] <= 1
 oneVehicle_constraint[1]: y[(6,_1)] + y[(6,_4)] <= 1
 oneVehicle_constraint[2]: y[(6,_2)] <= 1
 oneVehicle_constraint[3]: y[(6,_3)] <= 1
 R118: - x[(0,_0)] - x[(0,_1)] - x[(0,_2)] - x[(0,_3)] - x[(0,_4)]
   + sumPerDay[0] = 0
 constraint_deviationColSum[0]: deviationColSum[0] - sumPerDay[0] = -11.4
 R120: - x[(1,_0)] - x[(1,_1)] - x[(1,_2)] - x[(1,_3)] - x[(1,_4)]
   + sumPerDay[1] = 0
 constraint_deviationColSum[1]: deviationColSum[1] - sumPerDay[1] = -11.4
 R122: - x[(2,_0)] - x[(2,_1)] - x[(2,_2)] - x[(2,_3)] - x[(2,_4)]
   + sumPerDay[2] = 0
 constraint_deviationColSum[2]: deviationColSum[2] - sumPerDay[2] = -11.4
 R124: - x[(3,_0)] - x[(3,_1)] - x[(3,_2)] - x[(3,_3)] - x[(3,_4)]
   + sumPerDay[3] = 0
 constraint_deviationColSum[3]: deviationColSum[3] - sumPerDay[3] = -11.4
 R126: - x[(4,_0)] - x[(4,_1)] - x[(4,_2)] - x[(4,_3)] - x[(4,_4)]
   + sumPerDay[4] = 0
 constraint_deviationColSum[4]: deviationColSum[4] - sumPerDay[4] = -11.4
 R128: - x[(5,_0)] - x[(5,_1)] - x[(5,_2)] - x[(5,_3)] - x[(5,_4)]
   + sumPerDay[5] = 0
 constraint_deviationColSum[5]: deviationColSum[5] - sumPerDay[5] = -11.4
 R130: - x[(6,_0)] - x[(6,_1)] - x[(6,_2)] - x[(6,_3)] - x[(6,_4)]
   + sumPerDay[6] = 0
 constraint_deviationColSum[6]: deviationColSum[6] - sumPerDay[6] = -11.4
 constraint_deviation[(0,_0)]: - x[(0,_0)] + deviation[(0,_0)] = -4
 constraint_deviation[(1,_0)]: - x[(1,_0)] + deviation[(1,_0)] = -4
 constraint_deviation[(2,_0)]: - x[(2,_0)] + deviation[(2,_0)] = -4
 constraint_deviation[(3,_0)]: - x[(3,_0)] + deviation[(3,_0)] = -4
 constraint_deviation[(4,_0)]: - x[(4,_0)] + deviation[(4,_0)] = -4
 constraint_deviation[(5,_0)]: - x[(5,_0)] + deviation[(5,_0)] = -4
 constraint_deviation[(6,_0)]: - x[(6,_0)] + deviation[(6,_0)] = -4
 constraint_deviation[(0,_1)]: - x[(0,_1)] + deviation[(0,_1)] = -3.333
 constraint_deviation[(1,_1)]: - x[(1,_1)] + deviation[(1,_1)] = -3.333
 constraint_deviation[(2,_1)]: - x[(2,_1)] + deviation[(2,_1)] = -3.333
 constraint_deviation[(3,_1)]: - x[(3,_1)] + deviation[(3,_1)] = -3.333
 constraint_deviation[(4,_1)]: - x[(4,_1)] + deviation[(4,_1)] = -3.333
 constraint_deviation[(5,_1)]: - x[(5,_1)] + deviation[(5,_1)] = -3.333
 constraint_deviation[(6,_1)]: - x[(6,_1)] + deviation[(6,_1)] = -3.333
 constraint_deviation[(0,_2)]: - x[(0,_2)] + deviation[(0,_2)] = -4
 constraint_deviation[(1,_2)]: - x[(1,_2)] + deviation[(1,_2)] = -4
 constraint_deviation[(2,_2)]: - x[(2,_2)] + deviation[(2,_2)] = -4
 constraint_deviation[(3,_2)]: - x[(3,_2)] + deviation[(3,_2)] = -4
 constraint_deviation[(4,_2)]: - x[(4,_2)] + deviation[(4,_2)] = -4
 constraint_deviation[(5,_2)]: - x[(5,_2)] + deviation[(5,_2)] = -4
 constraint_deviation[(6,_2)]: - x[(6,_2)] + deviation[(6,_2)] = -4
 constraint_deviation[(0,_3)]: - x[(0,_3)] + deviation[(0,_3)] = -3.75
 constraint_deviation[(1,_3)]: - x[(1,_3)] + deviation[(1,_3)] = -3.75
 constraint_deviation[(2,_3)]: - x[(2,_3)] + deviation[(2,_3)] = -3.75
 constraint_deviation[(3,_3)]: - x[(3,_3)] + deviation[(3,_3)] = -3.75
 constraint_deviation[(4,_3)]: - x[(4,_3)] + deviation[(4,_3)] = -3.75
 constraint_deviation[(5,_3)]: - x[(5,_3)] + deviation[(5,_3)] = -3.75
 constraint_deviation[(6,_3)]: - x[(6,_3)] + deviation[(6,_3)] = -3.75
 constraint_deviation[(0,_4)]: - x[(0,_4)] + deviation[(0,_4)] = -4
 constraint_deviation[(1,_4)]: - x[(1,_4)] + deviation[(1,_4)] = -4
 constraint_deviation[(2,_4)]: - x[(2,_4)] + deviation[(2,_4)] = -4
 constraint_deviation[(3,_4)]: - x[(3,_4)] + deviation[(3,_4)] = -4
 constraint_deviation[(4,_4)]: - x[(4,_4)] + deviation[(4,_4)] = -4
 constraint_deviation[(5,_4)]: - x[(5,_4)] + deviation[(5,_4)] = -4
 constraint_deviation[(6,_4)]: - x[(6,_4)] + deviation[(6,_4)] = -4
Bounds
 deviationColSum[0] free
 deviationColSum[1] free
 deviationColSum[2] free
 deviationColSum[3] free
 deviationColSum[4] free
 deviationColSum[5] free
 deviationColSum[6] free
 deviation[(0,_0)] free
 deviation[(1,_0)] free
 deviation[(2,_0)] free
 deviation[(3,_0)] free
 deviation[(4,_0)] free
 deviation[(5,_0)] free
 deviation[(6,_0)] free
 deviation[(0,_1)] free
 deviation[(1,_1)] free
 deviation[(2,_1)] free
 deviation[(3,_1)] free
 deviation[(4,_1)] free
 deviation[(5,_1)] free
 deviation[(6,_1)] free
 deviation[(0,_2)] free
 deviation[(1,_2)] free
 deviation[(2,_2)] free
 deviation[(3,_2)] free
 deviation[(4,_2)] free
 deviation[(5,_2)] free
 deviation[(6,_2)] free
 deviation[(0,_3)] free
 deviation[(1,_3)] free
 deviation[(2,_3)] free
 deviation[(3,_3)] free
 deviation[(4,_3)] free
 deviation[(5,_3)] free
 deviation[(6,_3)] free
 deviation[(0,_4)] free
 deviation[(1,_4)] free
 deviation[(2,_4)] free
 deviation[(3,_4)] free
 deviation[(4,_4)] free
 deviation[(5,_4)] free
 deviation[(6,_4)] free
 Constant = 1
Binaries
 y[(0,_0)] y[(1,_0)] y[(2,_0)] y[(3,_0)] y[(4,_0)] y[(5,_0)] y[(6,_0)]
 y[(0,_1)] y[(1,_1)] y[(2,_1)] y[(3,_1)] y[(4,_1)] y[(5,_1)] y[(6,_1)]
 y[(0,_2)] y[(1,_2)] y[(2,_2)] y[(3,_2)] y[(4,_2)] y[(5,_2)] y[(6,_2)]
 y[(0,_3)] y[(1,_3)] y[(2,_3)] y[(3,_3)] y[(4,_3)] y[(5,_3)] y[(6,_3)]
 y[(0,_4)] y[(1,_4)] y[(2,_4)] y[(3,_4)] y[(4,_4)] y[(5,_4)] y[(6,_4)]
Generals
 x[(0,_0)] x[(1,_0)] x[(2,_0)] x[(3,_0)] x[(4,_0)] x[(5,_0)] x[(6,_0)]
 x[(0,_1)] x[(1,_1)] x[(2,_1)] x[(3,_1)] x[(4,_1)] x[(5,_1)] x[(6,_1)]
 x[(0,_2)] x[(1,_2)] x[(2,_2)] x[(3,_2)] x[(4,_2)] x[(5,_2)] x[(6,_2)]
 x[(0,_3)] x[(1,_3)] x[(2,_3)] x[(3,_3)] x[(4,_3)] x[(5,_3)] x[(6,_3)]
 x[(0,_4)] x[(1,_4)] x[(2,_4)] x[(3,_4)] x[(4,_4)] x[(5,_4)] x[(6,_4)]
 sumPerDay[0] sumPerDay[1] sumPerDay[2] sumPerDay[3] sumPerDay[4]
 sumPerDay[5] sumPerDay[6]
General Constraints
 constraint_absDeviationColSum[0]: absDeviationColSum[0] = ABS (
   deviationColSum[0] )
 constraint_absDeviationColSum[1]: absDeviationColSum[1] = ABS (
   deviationColSum[1] )
 constraint_absDeviationColSum[2]: absDeviationColSum[2] = ABS (
   deviationColSum[2] )
 constraint_absDeviationColSum[3]: absDeviationColSum[3] = ABS (
   deviationColSum[3] )
 constraint_absDeviationColSum[4]: absDeviationColSum[4] = ABS (
   deviationColSum[4] )
 constraint_absDeviationColSum[5]: absDeviationColSum[5] = ABS (
   deviationColSum[5] )
 constraint_absDeviationColSum[6]: absDeviationColSum[6] = ABS (
   deviationColSum[6] )
 constraint_absDeviation[(0,_0)]: absDeviation[(0,_0)] = ABS (
   deviation[(0,_0)] )
 constraint_absDeviation[(1,_0)]: absDeviation[(1,_0)] = ABS (
   deviation[(1,_0)] )
 constraint_absDeviation[(2,_0)]: absDeviation[(2,_0)] = ABS (
   deviation[(2,_0)] )
 constraint_absDeviation[(3,_0)]: absDeviation[(3,_0)] = ABS (
   deviation[(3,_0)] )
 constraint_absDeviation[(4,_0)]: absDeviation[(4,_0)] = ABS (
   deviation[(4,_0)] )
 constraint_absDeviation[(5,_0)]: absDeviation[(5,_0)] = ABS (
   deviation[(5,_0)] )
 constraint_absDeviation[(6,_0)]: absDeviation[(6,_0)] = ABS (
   deviation[(6,_0)] )
 constraint_absDeviation[(0,_1)]: absDeviation[(0,_1)] = ABS (
   deviation[(0,_1)] )
 constraint_absDeviation[(1,_1)]: absDeviation[(1,_1)] = ABS (
   deviation[(1,_1)] )
 constraint_absDeviation[(2,_1)]: absDeviation[(2,_1)] = ABS (
   deviation[(2,_1)] )
 constraint_absDeviation[(3,_1)]: absDeviation[(3,_1)] = ABS (
   deviation[(3,_1)] )
 constraint_absDeviation[(4,_1)]: absDeviation[(4,_1)] = ABS (
   deviation[(4,_1)] )
 constraint_absDeviation[(5,_1)]: absDeviation[(5,_1)] = ABS (
   deviation[(5,_1)] )
 constraint_absDeviation[(6,_1)]: absDeviation[(6,_1)] = ABS (
   deviation[(6,_1)] )
 constraint_absDeviation[(0,_2)]: absDeviation[(0,_2)] = ABS (
   deviation[(0,_2)] )
 constraint_absDeviation[(1,_2)]: absDeviation[(1,_2)] = ABS (
   deviation[(1,_2)] )
 constraint_absDeviation[(2,_2)]: absDeviation[(2,_2)] = ABS (
   deviation[(2,_2)] )
 constraint_absDeviation[(3,_2)]: absDeviation[(3,_2)] = ABS (
   deviation[(3,_2)] )
 constraint_absDeviation[(4,_2)]: absDeviation[(4,_2)] = ABS (
   deviation[(4,_2)] )
 constraint_absDeviation[(5,_2)]: absDeviation[(5,_2)] = ABS (
   deviation[(5,_2)] )
 constraint_absDeviation[(6,_2)]: absDeviation[(6,_2)] = ABS (
   deviation[(6,_2)] )
 constraint_absDeviation[(0,_3)]: absDeviation[(0,_3)] = ABS (
   deviation[(0,_3)] )
 constraint_absDeviation[(1,_3)]: absDeviation[(1,_3)] = ABS (
   deviation[(1,_3)] )
 constraint_absDeviation[(2,_3)]: absDeviation[(2,_3)] = ABS (
   deviation[(2,_3)] )
 constraint_absDeviation[(3,_3)]: absDeviation[(3,_3)] = ABS (
   deviation[(3,_3)] )
 constraint_absDeviation[(4,_3)]: absDeviation[(4,_3)] = ABS (
   deviation[(4,_3)] )
 constraint_absDeviation[(5,_3)]: absDeviation[(5,_3)] = ABS (
   deviation[(5,_3)] )
 constraint_absDeviation[(6,_3)]: absDeviation[(6,_3)] = ABS (
   deviation[(6,_3)] )
 constraint_absDeviation[(0,_4)]: absDeviation[(0,_4)] = ABS (
   deviation[(0,_4)] )
 constraint_absDeviation[(1,_4)]: absDeviation[(1,_4)] = ABS (
   deviation[(1,_4)] )
 constraint_absDeviation[(2,_4)]: absDeviation[(2,_4)] = ABS (
   deviation[(2,_4)] )
 constraint_absDeviation[(3,_4)]: absDeviation[(3,_4)] = ABS (
   deviation[(3,_4)] )
 constraint_absDeviation[(4,_4)]: absDeviation[(4,_4)] = ABS (
   deviation[(4,_4)] )
 constraint_absDeviation[(5,_4)]: absDeviation[(5,_4)] = ABS (
   deviation[(5,_4)] )
 constraint_absDeviation[(6,_4)]: absDeviation[(6,_4)] = ABS (
   deviation[(6,_4)] )
End
