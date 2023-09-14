# Import necessary packages
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import sys; sys.executable
import math
import csv
from datetime import datetime
import numpy as np




#### INPUT



### Gather relevant information of input files

    
## Gathering data from File 0: information on the hubs in each region

#Creating a dataframe hubInfo
hubInfo = pd.read_csv("data/File0.csv",header=0)

# Retrieve Hub names from dataframe and create a list that can be used at a later point
H = hubInfo['Hubs'].tolist()

# Number of hubs
numHubs = len(hubInfo)

# Last collection of previous week as list (per hub, same order); given in hubInfo File 0
lastCollection = hubInfo['LastCollection'].tolist()

# Next possible date of collection (no publVac/we); given in hubInfo File 0
nextCollection = hubInfo['NextCollection'].tolist()


## Gather information on prognosed amounts of collection: File 1

# Create dataframe with prognosis; columns: KW, rows: hubs
progVals = pd.read_csv("data/File1.csv",header=0)


## Find current week and necessary information in calendar: File 2 

# Create dataframe for values of File 2: calendar
cal = pd.read_csv("data/File2.csv",header=0)

# Find current date to later determine next Monday (start date of week to plan for)
currentDate = datetime.today().strftime('%Y-%m-%d')

# For comparison reasons: current date for output of the test version
currentDate = "2023-06-21"

# Find the index of the current date in File 2: calendar (as list of one value)
indexWeek = cal.index[cal['date'] == currentDate].tolist() 

# Get the index queried above as integer
indexWeek = indexWeek[0] 

# Initialise index to use as starting index for the week to plan
nextWeek = indexWeek 

# Make sure that the following week is used if current date is a monday (add 1 to produce a tuesday, see while-loop below)
if cal.at[indexWeek, 'wd'] == 1: 
    nextWeek = indexWeek + 1

# Find out next Monday in calendar to use the respective week as starting point in the calendar file (File 2)
while cal.at[nextWeek, 'wd'] != 1: 
    nextWeek += 1

# Find the selected date (Monday of week to prognose) to use it as key to later find the prognosed values that are used 
progDate = cal.at[nextWeek,'date'] 

# Number of days (given); Retrieve it from the dataframe created based on file1 after determining the prognosed date. The number of days 
# is always the last row of the file which is why H (the number of hubs) can be taken as index.
numDays = progVals.at[numHubs, progDate]

# Slice the calendar file to have only the values that are needed for the allocation of the next weeks 
# (given by the amount of days to prognose D): needed to find out where the next available days are
timeframeCal = cal[nextWeek:nextWeek+numDays] 


## Gather information of File 3: vehicle information

# Create a dataframe for vehicle information
vehicles = pd.read_csv("data/File3.csv",header=0) 

# Number of vehicles
numVehicles = len(vehicles) 

# Create a list of all vehicle IDs
vehicleIDs = vehicles['Vehicle-ID'].tolist()



### Calculate and store information that can be derived from given variables


## Determine information about collections

# Initialise a dictionary to store information on the collections per hub
Collections = {}

# Fill the dictionary with names of hubs an their respective number of collection days and collection tours

# Initialise control variable
j = 0

# Loop over hubs to store hub-specific information separately in the dictionary
for i in H:
    # Find vehicle that is associated with the respective hub
    vehicleID = hubInfo.at[j,'Vehicle-ID']  
    # Find out the index of the vehicle in dataframe "vehicles" to search for capacity later
    vehicleIDIndex = vehicles.index[vehicles['Vehicle-ID'] == vehicleID].tolist() 
    # Turn list value into integer value
    vehicleIDIndex = vehicleIDIndex[0] 
    # Calculate necessary amount of collection tours: round the value up to prevent having not enough collection tours 
    tours = math.ceil((progVals.at[j, progDate]*0.95)/vehicles.at[vehicleIDIndex,'Capacity'])
    # Calculate necessary amount of collection days: round the value up to prevent having not enough of collection days 
    days = math.ceil(tours/hubInfo.at[j,'Max'])
    # Find the average amount of collection tours per collection day and round the value to 3 digits; used for Objective 2: equal distribution of collection tours per hub
    average = round(tours/days,3)
    # Store variables in dictionary: the name of the current hub is used as key
    Collections[i] = {'Collection Tours': tours, 'Collection Days': days, 'Average Collection Tours': average}
    
    # Enable iteration over control variable
    j += 1



### Determine feasibility of model


## Find number and indices of days that are available for collection

# Initialise list to store days that allow for collection
availableDays = []

# Loop over days to find all possible collection days
for i in range(numDays):
    # Only append day as available day if it is not a day labelled as weekend or public vacation
    if (timeframeCal.at[nextWeek+i,'we'] == False) and (timeframeCal.at[nextWeek+i,'publVac'] == False):
        availableDays.append(i)

# Use the list of available days to store their amount as integer
numAvailDays = len(availableDays)


## Find hub groups: hubs that share the same vehicle

# Initialise dictionary to store information on each hub group
hubGroups = {}

# Initialise list to add each hub group as list (as many lists as there exist vehicles);
# will turn into list of lists
hubGroupsList = []

# Initialise control variable
l = 0

# Loop over the number of vehicles
for i in vehicleIDs:   
    # Initiate list to append on hub group list
    reusableList = []
    # Initiate entry in dictionary as list
    hubGroups[i] = {'Hub names':[]}
    # Initialise another control variable
    j = 0 
    
    # Iterate over number of hubs to see which hubs to add to the hub group
    for k in H:
        # If two or more hubs share a vehicle (that is, are assigned the same vehicle ID), 
        # then add them to the current list, when finished, continue with the next vehicle
        if hubInfo.at[j,'Vehicle-ID'] == i:
            # Append hub name to hub groups list 
            hubGroups[i]['Hub names'].append(k)
            # Append the index of the current hub to an empty list 
            # (in order to create a list that contains all indices of hubs that share the same vehicle)
            reusableList.append(j)
            
        # Enable iteration over control variable
        j += 1
        
    # Append the complete list of hub indices to the list of hub groups
    hubGroupsList.append([reusableList])
    # Enable iteration over control variable
    l += 1
    
        
## Find out whether there are enough collection days for one vehicle

# Initialise variable that states whether the model is feasible or not
feasibility = True

# Adding number of collection days as value for vehicles; loop over vehicles
for i in vehicleIDs:
    # Initialise a variable to sum up the collection days per vehicle    
    sumOfDays = 0
    # Initialise a variable to sum up the collection tours per vehicle
    sumOfTours = 0
    # Initialise a variable to store the average number of tours per collection day of the vehicle
    averageTours = 0
    
    # Loop over the number of hub groups
    for j in range(len(hubGroups[i]['Hub names'])):
        # Set a variable that is the name of the current hub
        hub = hubGroups[i]['Hub names'][j] 
        # Sum up the number of collection days
        sumOfDays += Collections[hub]['Collection Days']
        # Sum up the number of collection days
        sumOfTours += Collections[hub]['Collection Tours']
        
    # Store the calculated information: Sum of collection days  
    hubGroups[i].update({'Sum of collection days': sumOfDays })
    # Store the calculated information: Sum of collection tours
    hubGroups[i].update({'Sum of collection tours': sumOfTours })
    # Calculate the average number of tours per collection day of the vehicle
    averageTours = hubGroups[i]['Sum of collection tours']/hubGroups[i]['Sum of collection days']
    # Store the calculated information: Average collection tours
    hubGroups[i].update({'Average collection tours': averageTours })
    
    # If there are not enough collection days available for the specific vehicle, print that there have to be days added in the calendar file
    if sumOfDays > numAvailDays:
        print("Für diese Woche kann wegen Fahrzeug " + str(i) + " keine Lösung gefunden werden.") 
        # Set the feasibility variable to "False": Do not exit now but print all vehicles that are responsible for the infeasibility
        feasibility = False
        
# If one or more vehicles are short of available days, exit program: changing variables in File 0, File 2, or File 3 are necessary, 
# then the user should restart the program
if feasibility == False:
    sys.exit("Bitte entweder Sammeltage hinzufügen indem der Kalender bearbeitet wird, die maximale Anzahl der Touren pro Tag erhöhen oder ein Fahrzeug hinzufügen. Dann Programm neu starten.")




#### MODEL



### Set prerequisites for Objective 1: equal distribution throughout the week (deviation of sums)


## Sum up the number of collections in total to determine the mean

# Initialise variable that displays the number of collection tours        
totalCollections = 0

# Iterate over the number of hubs
for hub in H:
    # Per hub, add the amount of necessary collection tours to the amount of total collection tours
    totalCollections += Collections[hub]['Collection Tours']
       
# Calculate the mean: divide the sum of collection tours by the number of available days
meanTotal = totalCollections/numAvailDays



### Create prerequisite for Objective 3: equal distribution of spaces between collection days


## Creating penalty matrix for equal spaces between colections

# Initialise empty list; will turn into list of lists
penalties = []

# Iterate over number of hubs
for j in range (numHubs):
    # For each hub, append an empty list (create whole matrix): If there is a weekend or public holiday: state False, otherwise append a 10 (default value)
    penalties.append([])
    
    # Iterate over number of days to fill the matrix with initial values
    for i in range(numDays):
        penalties[j].append(10)

    
## Find the best possible collection days per hub (independent of hub groups)

# Loop over number of hubs
for j in range(numHubs):
    # Initiate variable that contains the name of the current hub
    hub = H[j]
    # Total sum of days that are of interest (excluding the last and the next possible collection)
    currentSum = nextCollection[j]-lastCollection[j]-1
    # Create an array that that contains all "indices" of the days of interest and 
    # reshape it such that these indices are the first row of the array
    table = np.arange((lastCollection[j]+1),nextCollection[j]).reshape((1,currentSum))
    # Initialise the second row: used to depict the availability of the respective day (if it allows for collection or not)
    secondRow = []
    
    # Loop over the number of days prior to the current week
    for i in range((abs(lastCollection[j]))-1):
        # Add zeros to all days: no collection possible before the start of the current week
        secondRow.append(0)
        
    # Loop over the number of days within the current week    
    for i in range(numDays):
        # If the day is not available for collection: 0; else use 1, indicates that collection is possible (np.array takes only variables of the same type such that these must be integer)
        if (timeframeCal.at[nextWeek+i,'we'] == True) or (timeframeCal.at[nextWeek+i,'publVac'] == True):
            secondRow.append(0)
        else:
            secondRow.append(1)
            
    # Append the second row to the existing array
    table = np.append(table, [secondRow], axis=0)
    # Create a third row that shall contain the collection days of the respective hub; start to fill it with zeros
    thirdRow = [0 for i in range(currentSum)]
    # Append the third row to the existing array
    table = np.append(table, [thirdRow], axis=0)    
    # Turn the array into a dataframe to enable easy handling in the later steps (could have been done from the start, theoretically)
    table = pd.DataFrame.from_records(table)
    # Set the first row that contains the indices of the three weeks combined as index axis to be able to access values by this index
    table.set_axis(table.iloc[0], axis=1, inplace=True)
    # Store the number of collections for the hub of concern in a new variable to be able to change its value without changing the original variable
    currentCollections = Collections[hub]['Collection Days']  
    
    # As long as there are collections to be assigned
    while currentCollections > 0:
        
        # Determine the current sum of days of concern
        currentSum = len(table.columns)
        # Find the first available day in the current week; use iloc. for "original index"
        firstAvailableDay = table.iloc[0, 0]
        
        # Loop over the number of current days of concern and stop as soon as the first available day is found
        for i in range(currentSum):
            # As soon as the day is indicated as possible collection day, set the variable and stop the for-loop
            if table.iloc[1, i] == 1:
                firstAvailableDay = table.iloc[0, i]
                break
        
        # Find the index of the last day where collection is allowed: iterate reversely to find the first day that is a 1 (collection allowed)
        lastAvailableDay = table.iloc[0, currentSum-1]

        # Reversely loop over the number of days of concern and count 
        for i in reversed(range(lastAvailableDay+1)):
            if table.at[1, i] == 1:
                break
            else:
                lastAvailableDay -= 1 
       
        # Find the ideal space between collections and round it
        idealSpace = round((currentSum-currentCollections)/(currentCollections+1))
        # Set the first collection day
        start = table.iloc[0,0] + idealSpace
        # Set the last collection of the current week to the same value to add the remaining spaces in the next step
        end = start
      
        # If there is no collection possible on the first ideal day
        if table.at[1, start] == 0:
            # Move the start day to the first day possible
            start = firstAvailableDay
        # Set the start day as collection day    
        table.at[2, start] = 1
        # One collection done: subtract it from the total number of collections in this hub
        currentCollections -= 1
    
        # Loop over the number of remaining collections (if only one collection, then it is already set)
        for i in range(currentCollections):
            # Add the ideal space to the variable "end" to obtain the ideal last day of collection
            end += idealSpace+1
    
        # If there are still collections left (if only one collection, then it is already set)
        if currentCollections > 0:

            # If the selected day is outside of the index range
            # This must be checked first to avoid an index access error!
            if end > (table.iloc[0, currentSum-1]):
                # Set the last collection day for this week to the last available day
                end = lastAvailableDay
            else:
                # If the selected day is not available
                if table.at[1, end] == 0:
                    # Set the last collection day for this week to the last available day
                    end = lastAvailableDay

            # Set the selected end date as collection day    
            table.at[2, end] = 1
            # Next collection done: subtract it from the total number of collections in this hub
            currentCollections -= 1
            
        # Transfer selected days into penalty matrix
        for i in table.iloc[0]:
            if table.at[2,i] == 1:
                penalties[j][i] = 0
        
        # Cut off the selected collection days to restart the process
        table = table.loc[:, start+1:end-1]


## Adjust penalty values in the penalty matrix    

# Loop over the number of hubs    
for j in range(numHubs):
    # Assign current hub name
    hub = H[j]
    # Find the index of the last collection day of the week
    lastCollectionDay = numDays-1
    # Set a control variable to find the last colleection day of the week
    n = numDays-1
    distanceToNextZero1 = 0
    distanceToNextZero2 = 0
    
    # Find the last collection day of the week by using the values in the penalty matrix
    while penalties[j][n] != 0:
        lastCollectionDay -= 1
        n -= 1
    
    # Loop over number of days
    for i in range(numDays):
        # Set a control variable that counts 
        l = i
        
        # If the current day is after the last collection, the distance to the next zero is at a maximum (take D to assure that it is always the highest)
        if i >= lastCollectionDay+1:
            distanceToNextZero1 = numDays
        
        # Else: Start to count the days until the next collection shall take place    
        else:
            # Set the counter to zero
            distanceToNextZero1 = 0  
            
            # As long as there is no zero found and the control variable is not out of range
            while penalties[j][l] != 0 and l<numDays-1:
                # Add one to the counter
                distanceToNextZero1 += 1
                l += 1
            
            # To avoid being out of the penalty matrix's range: add one more step if the while loop has not found the next zero yet
            if penalties[j][l] != 0:
                # Add one more to the counter
                distanceToNextZero1 += 1
                
        # Reset the control variable to start counting reversely
        l = i
        
        # If the current day is the first day in the week and is not a collection day: the reverse count must be at a maximum
        # (see above; use D to make sure it is always the highest value)
        if i == 0:
            # Set a new control variable for this purpose
            m = 0
            if penalties[j][i] != 0:
                while penalties[j][m] != 0 and m<numDays-1:
                    distanceToNextZero2 = numDays
                    m += 1
        # Else: proceed to start counting reversely
        else:
            # Set the counter to zero
            distanceToNextZero2 = 0 
            
            # As long as there is no zero found and the control variable is larger or equal to zero: add one to the counter
            while penalties[j][l] != 0 and l >= 0:
                distanceToNextZero2 += 1
                l -= 1
        
        # Transfer the smaller of both of the values to the penalty matrix
        penalties[j][i] = min(distanceToNextZero1, distanceToNextZero2)*10



### Create computational model

# Create model for gurobi
model = gp.Model()



### Add all variables the model for gurobi


## Create all decision variables

x = {}                 # Collection tours per day in each region: primary decision variable
y = {}                 # Binary variable indicating collection on each day in each region
deviation = {}         # Deviation per day per hub from average number of tours of the respective hub
absDeviation = {}      # Absolute value of the variable "deviation"
sumPerDay = {}         # Sum of tours in all hubs per day (column sum)
deviationColSum = {}    # Deviation of the total sum of tours per day (column sum) from average number of tours during the week
absDeviationColSum = {} # Absolute value of the deviationColSum


## Add all variables as objects of the 'gurobipy' module

# Loop over the number of hubs
for j in range(numHubs):
    # Loop over the number of days
    for i in range(numDays):
        # Adding the primary decision variable to model as integer (integer requirement according to Constraint 5 of the model)
        x[i, j] = model.addVar(vtype=GRB.INTEGER, name=f"x[{i, j}]") 
        # Adding an indicator variable that indicates whether there is a collection on the current day in the current hub or not: binary (requirement according to Constraint 6 of the model)
        # -> if 1: collection, if 0: no collection
        y[i, j] = model.addVar(vtype=GRB.BINARY, name=f"y[{i, j}]") 
        
        
## Add the total deviation variables to the model (Objective 1)

# Loop over number of days        
for i in range(numDays):
    # Add variable for deviation of the sum of colleections on this days from the mean of collections per day; lower bound is negative infinity
    deviationColSum[i] = model.addVar(vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY,  name=f"deviationColSum[{i}]")
    # Add variable that displays the absolute value of 'deviationColSum'
    absDeviationColSum[i] = model.addVar(vtype=GRB.CONTINUOUS, name=f"absDeviationColSum[{i}]")
    # Add variable that displays the sum of collection tours per day: is only needed for the output table
    sumPerDay[i] = model.addVar(vtype=GRB.INTEGER, name=f"sumPerDay[{i}]")

        
## Adding the deviation variables to the model (Objective 2):
  
# Loop over the number of hubs
for j in range(numHubs):
    # Loop over the number of days
    for i in range(numDays):          
        # As variables have a default lower bound of zero, the lower bound for the deviation variable has to be set to negative infinity
        deviation[i,j] = model.addVar(vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY,  name=f"deviation[{i,j}]")
        # Requirement above doesn't hold for absolute deviation: there, the default lower bound can be kept
        absDeviation[i,j] = model.addVar(vtype=GRB.CONTINUOUS, name=f"absDeviation[{i,j}]")
    


### Add constraints to the model 

   
## Add Constraints 1a, 1b and 1c: Possible values of x[i,j] must be considered; 
## by combination of constraint with y[i,j], Constraint 1c (0 as value is allowed) can be respected

# Loop over number of hubs
for j in range(numHubs):
    # Loop over number of days
    for i in range(numDays):
        # Adding Constraint 1a
        model.addConstr(x[i,j] >= hubInfo.at[j,'Min']  * y[i,j], f"min_constraint[{i,j}]")
        # Adding Constraint 1b
        model.addConstr(x[i,j] <= (hubInfo.at[j,'Max'])  * y[i,j], f"max_constraint[{i,j}]")        
        

## Add Constraints 2a and 2b: The number of collection tours and the number of collection days must be respected

# Loop over number of hubs
for j in range(numHubs):
    # Adding Constraint 2a
    model.addConstr(gp.quicksum(x[i,j] for i in range(numDays)) == Collections[hubInfo.at[j, 'Hubs']]['Collection Tours'], name=f"capacity_constraint_tours[{j}]")
    # Adding Constraint 2b
    model.addConstr(gp.quicksum(y[i,j] for i in range(numDays)) == Collections[hubInfo.at[j, 'Hubs']]['Collection Days'], name=f"capacity_constraint_days[{j}]")


## Add Constraint 3: Do not collect on days that are not allowed according to the calendar file

# Loop over number of days
for i in range(numDays):
    # Do not collect if either 'we' is True or if there is a public vacation day (i.e., 'publVac' is True)
    if (timeframeCal.at[nextWeek + i,'we'] == True) or (timeframeCal.at[nextWeek + i,'publVac'] == True):
        
        # Loop over number of hubs to add Constraint 3 to each hub
        for j in range(numHubs):
            model.addConstr((y[i, j] == 0), f"we/publVac_constraint[{i}]")


## Add Constraint 4: if hubs share the same vehicle, collection must not take place on the same day in these hubs

# Loop over number of hubs
#for j in range(H):
    # Loop over number of days
for i in range(numDays):
        # Loop over number of hub groups
    for k in range(len(hubGroupsList)):
        # Add the constraint such that there cannot be a collection if the vehicle is already used for another hub
        model.addConstr(((gp.quicksum(y[i, j] for j in hubGroupsList[k][0])) <= 1), f"oneVehicle_constraint[{k}]")
            

## Add Auxiliary Constraint 1: used to construct Objective 1

# Loop over number of days
for i in range(numDays):
    # Add Auxiliary Constraint 1a
    model.addConstr(sumPerDay[i] == gp.quicksum(x[i,j] for j in range(numHubs)))
    # Add Auxiliary Constraint 1b
    model.addConstr((deviationColSum[i] == (sumPerDay[i] - meanTotal)), f"constraint_deviationColSum[{i}]")    
    # Add Auxiliary Constraint 1c
    model.addGenConstrAbs(absDeviationColSum[i], deviationColSum[i], f"constraint_absDeviationColSum[{i}]")

## Add Auxiliary Constraint 2: used to construct Objective 2

# Loop over number of hubs
for j in range(numHubs):    
    # Loop over number of days
    for i in range(numDays):
        # Add Auxiliary Constraint 2a
        model.addConstr((deviation[i,j] == (x[i,j] - Collections[hubInfo.at[j, 'Hubs']]['Average Collection Tours'])),  f"constraint_deviation[{i,j}]")
        # Add Auxiliary Constraint 2b
        model.addGenConstrAbs(absDeviation[i,j], deviation[i,j], f"constraint_absDeviation[{i,j}]")



### Adding the objective function


## Adding each objective

# Initialise values of each objective
obj_1 = 0
obj_2 = 0
obj_3 = 0
  

## Objective 1: Equal distribution throughout the week (deviation of sums)

# Loop over number of days
for i in range(numDays):
    # Add the absolute value of total collection to the objective
    obj_1 += absDeviationColSum[i]

# Remove value of deviation for days where there is no collection
obj_1 -= (numDays-numAvailDays)*meanTotal 


## Objective 2: Equal distribution of tours per hub group

# Initialise variable that counts the days where no collection takes place
noCollection = 0

# Loop over number of hubs
for j in range(numHubs):
    # Loop over number of days
    for i in range(numDays):
        # Add the absolute value of each deviation to the objective 
        obj_2 += absDeviation[i,j]
   
    # Number of days that are not collected
    noCollection = numDays-Collections[hubInfo.at[j, 'Hubs']]['Collection Days']
    # Remove the deviation of days where there is no collection
    obj_2 -= noCollection*Collections[hubInfo.at[j, 'Hubs']]['Average Collection Tours']


## Objective 3: Equal distribution of collection days throughout the week

# Loop over number of hubs
for j in range(numHubs):
    # Loop over number of days
    for i in range(numDays):
        # Add the penalty value to the objective
        obj_3 += y[i,j] * penalties[j][i]


## Build the objective function by combining all three objectives
 
# Objective function as sum of all three objectives; Coefficients illustrate the weighting of each objective 
obj_fn = 0.8*obj_1 + 0.7*obj_2 + 0.4*obj_3

# Set the objective for the model
model.setObjective(obj_fn, GRB.MINIMIZE) 



### Set parameters for searching alternative feasible solutions

# Search intensely for alternative solutions and don't use the same solutions twice: pool search mode = 2
model.setParam(GRB.Param.PoolSearchMode, 2)

# Determines how many alternative solutions to retrieve: pool solutions = 3 + 1, because the first solution of the solution pool is always the optimal one
model.setParam(GRB.Param.PoolSolutions, 4)

# Optimise the model
model.optimize()




#### OUTPUT


## Save model

# LP-Format: for model browsing
model.write("model.lp")

# RLP-Format: same as LP-format, only different names for the variables (given by gurobi, not by program)
model.write("model.rlp")

# MPS-Format: to capture full model detail
model.write("model.mps")


## Retrieve the optimal as well as alternative solutions

# Only write output if optimal solution was found
if model.status == GRB.OPTIMAL:
    # Create a list that contains the headers for the output file
    headers = ["Hub", "Mo", "Di", "Mi", "Do", "Fr", "Sa", "So", "Touren pro Woche"]
    
    # Create a .csv file that is used to export the output of the program
    with open('data/Output.csv', 'w', newline='') as Output:
        # Assign the newly created file to the variable name "output"
        output = csv.writer(Output, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        # Write a row that contains the objective value of the optimal solution (for comparison reasons)
        output.writerow(["Zielwert:", round(model.objVal, 3)])
        # Write the headers as next row
        output.writerow(headers)
        # Add rows per hub that contain the hub name, the values per day as well as the sum of collections in the respective region
        for j in range(numHubs):
            outputList = [H[j]]
            for i in range(numDays):
                outputList.append(int(abs(x[i,j].x)))
            outputList.append(int(Collections[hubInfo.at[j, 'Hubs']]['Collection Tours']))    
            output.writerow(outputList)
        # Initialise a last row that contains the sum of collections per day    
        lastRow = ["Sum"]
        for i in range(numDays):
            lastRow.append(int(sumPerDay[i].x))
        lastRow.append(totalCollections)
        output.writerow(lastRow) 
        
        # Repeat procedure for writing down the optimal solution to also save the alternative solutions in the file
        for k in range(model.SolCount-1):
            model.setParam('SolutionNumber', k+1)
            output.writerow([" "])
            output.writerow(["Alternative", k+1])
            output.writerow(["Zielwert:", round(model.PoolObjVal, 3)])
            output.writerow(headers)
            for j in range(numHubs):
                outputList = [H[j]]
                for i in range(numDays):
                    outputList.append(int(abs(x[i,j].Xn)))
                outputList.append(int(Collections[hubInfo.at[j, 'Hubs']]['Collection Tours']))    
                output.writerow(outputList)
            lastRow = ["Summe"]
            for i in range(numDays):
                lastRow.append(int(sumPerDay[i].Xn))
            lastRow.append(totalCollections)
            output.writerow(lastRow)     
                
# If there is no feasible solution found, print it (should not happen due to exit point above (line 227)) 
else:
    print("Keine realisierbare Lösung gefunden.")