def ERxyz(Workspace, Runame, UMnamebase, Ufilemin, UFilesPerTimestep, TotalNumberOfNodes, Xmax, Ymax, SaveSpace):
    # This script calculates the Erosion Rate at each node.
    # A txt file is produced for each timestep listing ERs for all nodes in the same order as the xyz node files
    # A map of ERs is also prodcued for each timestep
    # The following should be changes before apply the script to a new model: (a) Workspace,
    # (b) Times, (c) zt1 location, (d) zt2 location, (e) UmapName, (f) df1 file name

    #Erosion rate is defined here for each node as:
    #The mean Uplift experienced by the node minus (the change in elevation divided by the total time between the adjacent timesteps)
    #There should be no negative values of ER in a detachment-limited model

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    # 2Import the data. This is the csv xyz file (from excel conversion after CHILD_to_xyz)
    Times = open(Workspace + Runame + 'timeslices.txt', "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    n = 0
    ntot = len(Timeslices)
    Ufilemax = Ufilemin + UFilesPerTimestep
    Ufilemin = Ufilemin

    while n < (ntot-1):
        Zt1time = Timeslices[n]
        Zt1location = Workspace + Runame + '_time' + str(Zt1time) + '.csv'
        Zt1data = pd.read_csv(Zt1location)
        Zt1 = Zt1data['z']

        Zt2time = Timeslices[n+1]
        Zt2location = Workspace + Runame + '_time' + str(Zt2time) + '.csv'
        Zt2data = pd.read_csv(Zt2location)
        Zt2 = Zt2data['z']

        DeltaZ = []
        DeltaZ = Zt2 - Zt1

        Ufileends = list(range(Ufilemin, Ufilemax))

        df = pd.DataFrame(index=np.arange(TotalNumberOfNodes))

        for f in Ufileends:
            f3 = str(f).zfill(3)
            UmapName = UMnamebase + str(f3)
            Umap = Workspace + UmapName
            Udata = pd.read_csv(Umap, header=None, names=['A'])
            us = Udata['A']
            df[f3] = us
            print('Uplift Map' + f3 +'Added to Mean Calculation')

        MeanU = df.mean(axis=1)
        print ('Mean U Calculated for Timesteps' + str(n))

        Time = int(Timeslices[n + 1]) - int(Timeslices[n])

        ERs = []
        ERs = MeanU - (DeltaZ/Time)

        # Add ERs to PD df with Z xy nodes
        Zt1data['ER'] = ERs
        ERdataframeOutputName = Workspace + 'ERxyz' + str(Timeslices[n]) + '.csv'
        Zt1data.to_csv(ERdataframeOutputName)

        # Save ERs to txt file
        OutERFileName = Workspace + 'ER' + str(Timeslices[n]) + '.txt'
        with open(OutERFileName, 'w') as f:
            for line in ERs:
                f.writelines(str(line))
                f.writelines('\n')

        # ERFileName = open(OutERFileName, 'w')
        # ERFileName.write(str(ERs))
        # ERFileName.close()

        n = n + 1

        Ufilemin = Ufilemax + 1
        Ufilemax = Ufilemax + UFilesPerTimestep

    for i in range(0, len(Timeslices) - 1):  # If we want to do that for each time slice, and plot one after the other
        # for i in range (1, 2):               # If we just want to focus on the second time slice
        time = Timeslices[i]
        EroRates = pd.read_csv(Workspace + "ERxyz" + str(time) + ".csv")
        ER = EroRates['ER']
        print("Processing time slice for time =" + str(time) + " years.")
        #df1 = pd.read_csv(Workspace + Runame + "_time" + str(time) + '.csv')
        x = EroRates['x']
        y = EroRates['y']

        # Map
        fig1 = plt.figure(str(i) + 'a')
        plt.scatter(x, y, c=(ER*1000))
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title('time = ' + str(time) + 'yr')

        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)

        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
        plt.axis([0, Xmax, 0, Ymax])
        cbar.set_label('Erosion Rate (mm/yr)')

        # Write the figure to file
        plt.savefig(SaveSpace + Runame + str(i) + '.png')

def ERcs(Workspace, Runame, PlotXmin, PlotXmax, PlotYmin, PlotYmax, SaveSpace):
    # The following should be changed before applying the file to a new model:
    # (1) workspace, (2) Times, (3)

    # 1Import relevant modules, check extensions available, and create environment workspace
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    # 2Import the data. This is the txt file of timeslices from the CHILD_to_xyz script
    Times = open(Workspace + Runame + 'timeslices2run.txt', "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    # 3Import the data. The .csv ER data, and the .csv MDD location data
    ntotal = len(Timeslices)
    n = 0

    while n < ntotal -1:
        # Take the ER file and the MDD file for each timestep and make into a point file (.shp)
        ERxyzfilename = 'ERxyz' + Timeslices[n] + '.csv'
        MDDfilename = 'MDDsa' + Timeslices[n] + '.csv'

        # Convert both to a pd dataframe
        dfNodes = pd.read_csv(Workspace + ERxyzfilename)
        dfMDD = pd.read_csv(Workspace + MDDfilename)

        Nodecoordinates = dfNodes[['x', 'y']].values.tolist()
        MDDcoordinates = dfMDD[['POINT_X', 'MEAN_POINT_Y']].values.tolist()
        NodeYcoords = dfNodes['y'].values.tolist()
        MDDYcoords = dfMDD['MEAN_POINT_Y'].values.tolist()

        Nodecoordinates = np.array(Nodecoordinates)
        MDDcoordinates = np.array(MDDcoordinates)

        NodeMDDdistances = []
        NodeMDDdirection = []

        # Here the distance between each node and each point of the MDD is calculated. The smallest distance is taken as the distance between the node and the MDD.
        for i in Nodecoordinates:
            distances = np.linalg.norm(MDDcoordinates - i, axis=1)
            min_index = np.argmin(distances)
            min_dist = distances[min_index]
            direction = MDDYcoords[min_index] - i[1]
            NodeMDDdistances.append(min_dist)
            NodeMDDdirection.append(direction)

        # Add min distances to PD df with ERs and xyz nodes
        dfNodes['DistanceMDD'] = NodeMDDdistances
        dfNodes['DirectionMDD'] = NodeMDDdirection
        ERMDDdataframeOutputName = 'ERMDDxyz' + str(Timeslices[n]) + '.csv'
        dfNodes.to_csv(Workspace + ERMDDdataframeOutputName)

        # Save ERs to txt file
        OutNodeDistMDD = Workspace + 'NodeDistMDD' + str(Timeslices[n]) + '.txt'
        with open(OutNodeDistMDD, 'w') as f:
            for line in NodeMDDdistances:
                f.writelines(str(line))
                f.writelines('\n')

        print(str(Timeslices[n]))
        n = n + 1

    n = 0
    ntotal = len(Timeslices)
    while n < ntotal -1:
        ERdataFile = Workspace + 'ERMDDxyz' + str(Timeslices[n]) + '.csv'
        ERdata = pd.read_csv(ERdataFile)
        # ERdata['NEAR_ANGLE'].replace(ERdata['NEAR_ANGLE']<0, -1)
        ERdata['DirectionMDD'] = ERdata['DirectionMDD'].mask(ERdata['DirectionMDD'] < 0, -1)
        ERdata['DirectionMDD'] = ERdata['DirectionMDD'].mask(ERdata['DirectionMDD'] >= 0, 1)
        ERdata['MDDdist'] = ERdata['DistanceMDD'] * ERdata['DirectionMDD']
        ERdata['MDDdist01km'] = (ERdata['MDDdist'] / 1000).round(1)

        ERdistMeans = ERdata.groupby('MDDdist01km')['ER'].mean()
        ERdistances = ERdata.MDDdist01km.unique()
        ERdistances.sort()
        # stats.binned_statistic(Bins, ERdata[''])

        Dist_km = ERdata['MDDdist'] / 1000

        # Plot ER Data
        fig1, (ax1) = plt.subplots(nrows=1, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1]})
        ax1.scatter(Dist_km, (ERdata['ER']*1000), color='darkgrey', s=0.1, zorder=2)
        ax1.plot(ERdistances, (ERdistMeans*1000), color='k', zorder=3)
        # ax1.plot(Dist_kmbZ, binnedZ['MIN_grid_code'], color='b')
        # ax1.plot(Dist_kmbZ, binnedZ['MAX_grid_code'], color='r')
        ax1.axis([PlotXmin, PlotXmax, PlotYmin, PlotYmax])
        ax1.set_ylabel('Erosion Rate (mm/yr)')
        ax1.set_xlabel('[Distal]           Distance from MDD (km)           [Proximal]')
        ax1.axvline(x=0, color='lightgrey', zorder=1)
        ax1.axhline(y=0, color='lightgrey', zorder=1)
        fig1.suptitle(str(Timeslices[n]))

        # Write the figure to file
        plt.savefig(SaveSpace + Runame + str(Timeslices[n]) + '.png')

        print(Timeslices[n])
        n = n + 1

def ERxyzUG(Workspace, Runame, Umap):
    # This script calculates the Erosion Rate at each node.
    # A txt file is produced for each timestep listing ERs for all nodes in the same order as the xyz node files
    # A map of ERs is also prodcued for each timestep
    # The following should be changes before apply the script to a new model: (a) Workspace,
    # (b) Times, (c) zt1 location, (d) zt2 location, (e) UmapName, (f) df1 file name

    #This script is for when there is no advection (i.e. when the umap remains the same with time)

    #Erosion rate is defined here for each node as:
    #The mean Uplift experienced by the node minus (the change in elevation divided by the total time between the adjacent timesteps)
    #There should be no negative values of ER in a detachment-limited model

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    # 2Import the data. This is the csv xyz file (from excel conversion after CHILD_to_xyz)
    Times = open(Workspace + Runame + 'timeslices.txt', "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    n = 0
    ntot = len(Timeslices)


    while n < (ntot-1):
        Zt1time = Timeslices[n]
        Zt1location = Workspace + Runame + '_time' + str(Zt1time) + '.csv'
        Zt1data = pd.read_csv(Zt1location)
        Zt1 = Zt1data['z']

        Zt2time = Timeslices[n+1]
        Zt2location = Workspace + Runame + '_time' + str(Zt2time) + '.csv'
        Zt2data = pd.read_csv(Zt2location)
        Zt2 = Zt2data['z']

        DeltaZ = []
        DeltaZ = Zt2 - Zt1

        #Here we preform the calculation on DeltaZ and Ufile
        Udata = pd.read_csv(Umap, header=None, names=['A'])
        us = Udata['A']
        print('Uplift Map Added to Mean Calculation')
        print ('Mean U Calculated for Timesteps' + str(n))

        Time = int(Timeslices[n + 1]) - int(Timeslices[n])

        ERs = []
        ERs = us - (DeltaZ/Time)

        # Add ERs to PD df with Z xy nodes
        Zt1data['ER'] = ERs
        ERdataframeOutputName = Workspace + 'ERxyz' + str(Timeslices[n]) + '.csv'
        Zt1data.to_csv(ERdataframeOutputName)

        # Save ERs to txt file
        OutERFileName = Workspace + 'ER' + str(Timeslices[n]) + '.txt'
        with open(OutERFileName, 'w') as f:
            for line in ERs:
                f.writelines(str(line))
                f.writelines('\n')

        # ERFileName = open(OutERFileName, 'w')
        # ERFileName.write(str(ERs))
        # ERFileName.close()

        n = n + 1

def ERxyzOneUmap(Workspace, Runame, Ufilemin, UFilesPerTimestep, TotalNumberOfNodes, Xmax, Ymax, SaveSpace):
    # This script calculates the Erosion Rate at each node.
    # A txt file is produced for each timestep listing ERs for all nodes in the same order as the xyz node files
    # A map of ERs is also prodcued for each timestep
    # The following should be changes before apply the script to a new model: (a) Workspace,
    # (b) Times, (c) zt1 location, (d) zt2 location, (e) UmapName, (f) df1 file name

    #Erosion rate is defined here for each node as:
    #The mean Uplift experienced by the node minus (the change in elevation divided by the total time between the adjacent timesteps)
    #There should be no negative values of ER in a detachment-limited model

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    # 2Import the data. This is the csv xyz file (from excel conversion after CHILD_to_xyz)
    Times = open(Workspace + Runame + 'timeslices.txt', "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    n = 0
    ntot = len(Timeslices)
    Ufilemax = Ufilemin + UFilesPerTimestep
    Ufilemin = Ufilemin

    while n < (ntot - 1):
        Zt1time = Timeslices[n]
        Zt1location = Workspace + Runame + '_time' + str(Zt1time) + '.csv'
        Zt1data = pd.read_csv(Zt1location)
        Zt1 = Zt1data['z']

        Zt2time = Timeslices[n + 1]
        Zt2location = Workspace + Runame + '_time' + str(Zt2time) + '.csv'
        Zt2data = pd.read_csv(Zt2location)
        Zt2 = Zt2data['z']

        DeltaZ = []
        DeltaZ = Zt2 - Zt1

        Ufileends = list([1])

        df = pd.DataFrame(index=np.arange(TotalNumberOfNodes))

        for f in Ufileends:
            f3 = str(f).zfill(3)
            UmapName = 'UM' + Runame + str(f3)
            Umap = Workspace + UmapName
            Udata = pd.read_csv(Umap, header=None, names=['A'])
            us = Udata['A']
            df[f3] = us
            print('Uplift Map' + f3 + 'Added to Mean Calculation')

        MeanU = df.mean(axis=1)
        print ('Mean U Calculated for Timesteps' + str(n))

        Time = int(Timeslices[n + 1]) - int(Timeslices[n])

        ERs = []
        ERs = MeanU - (DeltaZ / Time)

        # Add ERs to PD df with Z xy nodes
        Zt1data['ER'] = ERs
        ERdataframeOutputName = Workspace + 'ERxyz' + str(Timeslices[n]) + '.csv'
        Zt1data.to_csv(ERdataframeOutputName)

        # Save ERs to txt file
        OutERFileName = Workspace + 'ER' + str(Timeslices[n]) + '.txt'
        with open(OutERFileName, 'w') as f:
            for line in ERs:
                f.writelines(str(line))
                f.writelines('\n')

        # ERFileName = open(OutERFileName, 'w')
        # ERFileName.write(str(ERs))
        # ERFileName.close()

        n = n + 1

    for i in range(0, len(Timeslices) - 1):  # If we want to do that for each time slice, and plot one after the other
        # for i in range (1, 2):               # If we just want to focus on the second time slice
        time = Timeslices[i]
        EroRates = pd.read_csv(Workspace + "ERxyz" + str(time) + ".csv")
        ER = EroRates['ER']
        print("Processing time slice for time =" + str(time) + " years.")
        # df1 = pd.read_csv(Workspace + Runame + "_time" + str(time) + '.csv')
        x = EroRates['x']
        y = EroRates['y']

        # Map
        fig1 = plt.figure(str(i) + 'a')
        plt.scatter(x, y, c=(ER * 1000))
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title('time = ' + str(time) + 'yr')

        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)

        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
        plt.axis([0, Xmax, 0, Ymax])
        cbar.set_label('Erosion Rate (mm/yr)')

        # Write the figure to file
        plt.savefig(SaveSpace + Runame + str(i) + '.png')


def ERxyzNoUmap(Workspace, Runame, UniformUrate, Xmax, Ymax, SaveSpace):
    # This script calculates the Erosion Rate at each node.
    # A txt file is produced for each timestep listing ERs for all nodes in the same order as the xyz node files
    # A map of ERs is also prodcued for each timestep
    # The following should be changes before apply the script to a new model: (a) Workspace,
    # (b) Times, (c) zt1 location, (d) zt2 location, (e) UmapName, (f) df1 file name

    #Erosion rate is defined here for each node as:
    #The mean Uplift experienced by the node minus (the change in elevation divided by the total time between the adjacent timesteps)
    #There should be no negative values of ER in a detachment-limited model

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    # 2Import the data. This is the csv xyz file (from excel conversion after CHILD_to_xyz)
    Times = open(Workspace + Runame + 'timeslices.txt', "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    n = 0
    ntot = len(Timeslices)

    while n < (ntot - 1):
        Zt1time = Timeslices[n]
        Zt1location = Workspace + Runame + '_time' + str(Zt1time) + '.csv'
        Zt1data = pd.read_csv(Zt1location)
        Zt1 = Zt1data['z']

        Zt2time = Timeslices[n + 1]
        Zt2location = Workspace + Runame + '_time' + str(Zt2time) + '.csv'
        Zt2data = pd.read_csv(Zt2location)
        Zt2 = Zt2data['z']

        DeltaZ = []
        DeltaZ = Zt2 - Zt1

        MeanU = UniformUrate
        print ('U set for' + str(n))

        Time = int(Timeslices[n + 1]) - int(Timeslices[n])

        ERs = []
        ERs = MeanU - (DeltaZ / Time)

        # Add ERs to PD df with Z xy nodes
        Zt1data['ER'] = ERs
        ERdataframeOutputName = Workspace + 'ERxyz' + str(Timeslices[n]) + '.csv'
        Zt1data.to_csv(ERdataframeOutputName)

        # Save ERs to txt file
        OutERFileName = Workspace + 'ER' + str(Timeslices[n]) + '.txt'
        with open(OutERFileName, 'w') as f:
            for line in ERs:
                f.writelines(str(line))
                f.writelines('\n')

        # ERFileName = open(OutERFileName, 'w')
        # ERFileName.write(str(ERs))
        # ERFileName.close()

        n = n + 1

    for i in range(0, len(Timeslices) - 1):  # If we want to do that for each time slice, and plot one after the other
        # for i in range (1, 2):               # If we just want to focus on the second time slice
        time = Timeslices[i]
        EroRates = pd.read_csv(Workspace + "ERxyz" + str(time) + ".csv")
        ER = EroRates['ER']
        print("Processing time slice for time =" + str(time) + " years.")
        # df1 = pd.read_csv(Workspace + Runame + "_time" + str(time) + '.csv')
        x = EroRates['x']
        y = EroRates['y']

        # Map
        fig1 = plt.figure(str(i) + 'a')
        plt.scatter(x, y, c=(ER * 1000))
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title('time = ' + str(time) + 'yr')

        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)

        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
        plt.axis([0, Xmax, 0, Ymax])
        cbar.set_label('Erosion Rate (mm/yr)')

        # Write the figure to file
        plt.savefig(SaveSpace + Runame + str(i) + '.png')