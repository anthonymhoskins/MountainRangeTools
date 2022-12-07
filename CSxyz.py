#For plotting CS of entire model
def CS(Workspace, Runame, Ymax):

    #import the required modeules
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # 2Import the data. This is the txt file of timeslices from the CHILD_to_xyz script.
    Times = open(Workspace + Runame + "timeslicesCS.txt", "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    ntotal = len(Timeslices)
    n = 0

    fig1, (ax1) = plt.subplots(nrows=1, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1]})

    while n < ntotal:
        xyz = pd.read_csv(Workspace + Runame + '_time' + Timeslices[n] + '.csv')
        xyz = xyz.round({'y':-2})
        groups = xyz.groupby('y')['z'].mean()
        z = list(groups)
        y = groups.index.tolist()
        ax1.plot(y, z, color='k', alpha=((1/ntotal) * (n+1)))
        ax1.axis([0, Ymax, 0, 500])
        ax1.set_ylabel('Elevation (m)')
        ax1.set_xlabel('Cross Section')
        ax1.set_title(Runame)
        print(n)
        n=n+1

    #Write the figure to file
    plt.savefig(Workspace + Runame + 'CS.png')

#For plotting CS across active area of uplift only when advection is induced
def CSstandardFaultAD(Workspace, Runame, WW, ADrate, NonUpliftY, SaveName):

    #import the required modeules
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # 2Import the data. This is the txt file of timeslices from the CHILD_to_xyz script.
    Times = open(Workspace + Runame + "timeslicesCS.txt", "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    ntotal = len(Timeslices)
    n = 0

    fig1, (ax1) = plt.subplots(nrows=1, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1]})

    while n < ntotal:
        xyz = pd.read_csv(Workspace + Runame + '_time' + Timeslices[n] + '.csv')
        xyz = xyz.round({'y':-2})
        groups = xyz.groupby('y')['z'].mean()
        z = list(groups)
        y = groups.index.tolist()
        for i in range(len(y)):
            y[i] -= NonUpliftY - (ADrate * ((int(Timeslices[n]))-200000000))
        ax1.plot(y, z, color='k', alpha=((1/ntotal) * (n+1)))
        ax1.axis([0, WW, 0, 500])
        ax1.set_ylabel('Elevation (m)')
        ax1.set_xlabel('Cross Section (m)')
        ax1.set_title('Cross Section' + ' ' + Runame)
        print(n)
        n=n+1

    #Write the figure to file
    plt.savefig(Workspace + Runame + SaveName + '.png')
