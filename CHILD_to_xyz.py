

def CHILDxyz(CHILDlocation,Runame,SaveLocation):
    import pandas as pd  # We will use pandas later
    import os
    import csv

    runname = Runame # Specify the name of your run

    fz = open(CHILDlocation + runname + ".z", 'r')  # open .z file
    linesz = fz.readlines()  # read in the data
    n_lines = len(linesz)  # get the number of lines

    fxy = open(CHILDlocation + runname + ".nodes", 'r')  # open the .nodes file
    linesxy = fxy.readlines()  # read in the data

    # Initialisation of variables for the loop
    keepgoing = 1
    rowref = 0  # This is the row reference number, starting at 0. The number of node will be added at each time step
    # to start the procedure of extracting the relevant data from that point. This will also be used to run the test
    # to stop the keepgoing loop: if rowref becomes > to the number of lines in the file --> stop (keepgoing = 0)
    timeslices = [];  # We will record the different time slices in this array, to be able to use them for plotting later

    while keepgoing:
        x = [];
        y = [];
        z = [];
        b = [];
        data = [];
        df = [];  # Initialisation of arrays / variables
        d = (linesz[rowref])
        ds = round(float(d))
        time = int(ds)  # Read the time from the row = rowref in the z array
        print("Processing time slice for time = " + str(time) + " years.")
        nn = int(
            linesz[rowref + 1])  # Read the number of nodes for this time slice from the row = rowref+1 in the z array
        for i in range(rowref + 2, rowref + nn + 2):  # Get the data for all the nodes in this time slice
            linesxy_array = linesxy[i].strip().split(
                " ")  # Split the lines from the .nodes file into elements using the space as the delimiter
            x.append(float(linesxy_array[0]))  # x is the first element
            y.append(float(linesxy_array[1]))  # y is the second element
            b.append(int(linesxy_array[3]))  # b is the fourth element
            z.append(float(linesz[i]))  # take z from the linesz array

        timeslices.append(int(time))  # Save the different times here for plotting later

        data = {'x': x, 'y': y, 'z': z, 'b': b}  # put all the data for this time slice in an array
        df = pd.DataFrame(data)  # convert the dictionary into a data frame with pandas
        # df # This is our xyzb data!

        outfilenamea = SaveLocation + runname + "_time" + str(
            time) + ".txt"  # Writing the data in a csv file that has a name made of run name + time considered
        outfilea = open(outfilenamea, 'w', )
        outfilea.write(df.to_string())
        outfilea.close()

        outfilenameb = SaveLocation + runname + "_time" + str(time) + ".csv"
        df.to_csv(outfilenameb)

        # Done! Let's move to the next time step!
        rowref = rowref + nn + 2
        if rowref >= n_lines:
            keepgoing = 0
            print("Job done!")

    #timeslices = timeslices.replace["[", ""]
    timeslices = str(timeslices)[1:-1]
    outfilenamec = SaveLocation + runname + "timeslices.txt"
    outfilec = open(outfilenamec, 'w')
    outfilec.write(str(timeslices))
    outfilec.close()

def Mapxyz(CHILDlocation, Runame, Ymax, SaveSpace):
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    import time
    import pandas as pd

    Timefile = open(CHILDlocation + Runame + 'timeslices.txt', "r")
    content = Timefile.read()
    timeslices = content.split(", ")
    Timefile.close()

    for i in range(0, len(timeslices)):  # If we want to do that for each time slice, and plot one after the other
        # for i in range (1, 2):               # If we just want to focus on the second time slice
        df2 = []
        time = timeslices[i]
        print("Processing time slice for time =" + str(time) + " years.")
        df2 = pd.read_csv(CHILDlocation + Runame + "_time" + str(time) + ".txt", sep="\s+")
        x = df2.x
        y = df2.y
        z = df2.z

        # Map
        fig1 = plt.figure(str(i) + 'a')
        plt.scatter(x, y, c=z)
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.title('time = ' + str(time) + 'yr')
        cbar = plt.colorbar()
        plt.axis([0, 50000, 0, Ymax])
        cbar.set_label('Elevation (m)')

        # Write the figure to file
        plt.savefig(SaveSpace + Runame + str(i) + '.png')

def MapxyzAD(CHILDlocation, Runame, Ymax, SaveSpace, AD):
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    import time
    import pandas as pd

    Timefile = open(CHILDlocation + Runame + 'timeslices.txt', "r")
    content = Timefile.read()
    timeslices = content.split(", ")
    Timefile.close()

    for i in range(0, len(timeslices)):  # If we want to do that for each time slice, and plot one after the other
        # for i in range (1, 2):               # If we just want to focus on the second time slice
        df2 = []
        time = timeslices[i]
        print("Processing time slice for time =" + str(time) + " years.")
        df2 = pd.read_csv(CHILDlocation + Runame + "_time" + str(time) + ".txt", sep="\s+")
        x = (df2.x)/1000
        y = (df2.y)/1000
        z = df2.z

        YminNow = 0 - ((AD/1000)*(int(time)-200000000))
        YmaxNow = Ymax - ((AD/1000)*(int(time)-200000000))
        Yytext = 19 - (AD/1000)*(int(time)-200000000)

        # Map
        fig1 = plt.figure(str(i) + 'a')
        fig1.figsize=(8,6)
        plt.scatter(x, y, c=z, cmap='viridis', vmin=0, vmax=600)
        plt.xlabel('X (km)')
        #plt.ylabel('Y (km)')
        plt.text(-6, Yytext, 'Y (km)', rotation='vertical')
        plt.title('time = ' + str(time) + 'yr')
        cbar = plt.colorbar()
        plt.axis([0, 50, YminNow, YmaxNow])
        cbar.set_label('Elevation (m)')
        # Write the figure to file
        plt.savefig(SaveSpace + Runame + str(i) + '.png')

def MapxyzNOAD(CHILDlocation, Runame, Ymax, SaveSpace):
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    import time
    import pandas as pd

    Timefile = open(CHILDlocation + Runame + 'timeslices.txt', "r")
    content = Timefile.read()
    timeslices = content.split(", ")
    Timefile.close()

    for i in range(0, len(timeslices)):  # If we want to do that for each time slice, and plot one after the other
        # for i in range (1, 2):               # If we just want to focus on the second time slice
        df2 = []
        time = timeslices[i]
        print("Processing time slice for time =" + str(time) + " years.")
        df2 = pd.read_csv(CHILDlocation + Runame + "_time" + str(time) + ".txt", sep="\s+")
        x = (df2.x)/1000
        y = (df2.y)/1000
        z = df2.z

        YminNow = 0
        YmaxNow = Ymax

        # Map
        fig1 = plt.figure(str(i) + 'a')
        fig1.figsize=(8,6)
        plt.scatter(x, y, c=z, cmap='viridis', vmin=0, vmax=600)
        plt.xlabel('X (km)')
        #plt.ylabel('Y (km)')
        plt.text(-6, 19, 'Y (km)', rotation='vertical')
        plt.title('time = ' + str(time) + 'yr')
        cbar = plt.colorbar()
        plt.axis([0, 50, YminNow, YmaxNow])
        cbar.set_label('Elevation (m)')
        # Write the figure to file
        plt.savefig(SaveSpace + Runame + str(i) + '.png')
