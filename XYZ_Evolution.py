def Zevolution(BLOCKWorkspace, BLOCKRuname, UGWorkspace, UGRuname, UGADWorkspace, UGADRuname, StormInterval, InitialFaultPositionBLOCK, InitialFaultPositionUG, InitialFaultPositionUGAD, ProprateBLOCK, ProprateUG, ProprateUGAD, WedgeWidth, BLOCKTimeAdjust, UGTimeAdjust, UGADTimeAdjust, Figure, FigureTitle):

    #Import the dependent modules
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    #Load in the .vols files for the plot - this is the landscape volume recorded for each storm
    MassBLOCK = pd.read_csv(BLOCKWorkspace + BLOCKRuname + '.vols')
    MassUG = pd.read_csv(UGWorkspace + UGRuname + '.vols')
    MassUGAD = pd.read_csv(UGADWorkspace + UGADRuname + '.vols')

    #Find the number of entries in the files
    MassBLOCKLen = len(MassBLOCK)
    MassUGLen = len(MassUG)
    MassUGADLen = len(MassUGAD)

    #Create a list giving the range of times when landscape volume is recorded
    TimeRangeBLOCK = list(range(0,MassBLOCKLen))
    TimeRangeUG = list(range(MassBLOCKLen, (MassBLOCKLen+MassUGLen)))
    TimeRangeUGAD = list(range((MassBLOCKLen+MassUGLen), (MassBLOCKLen+MassUGLen+MassUGADLen)))
    ST_STDUR = StormInterval
    TimeRangeBLOCKArray = np.array(TimeRangeBLOCK)
    TimeRangeUGArray = np.array(TimeRangeUG)
    TimeRangeUGADArray = np.array(TimeRangeUGAD)
    TimeBLOCK = (TimeRangeBLOCKArray * ST_STDUR)
    TimeUG = (TimeRangeUGArray * ST_STDUR)
    TimeUGAD = (TimeRangeUGADArray * ST_STDUR)

    #Back to the future Z
    #Load the xyz data from csv form
    TimesBLOCK = open(BLOCKWorkspace + BLOCKRuname + 'timeslices.txt', "r")
    ContentBLOCK = TimesBLOCK.read()
    TimeslicesBLOCK = ContentBLOCK.split(", ")
    TimesBLOCK.close()
    print(TimeslicesBLOCK)
    TimesUG = open(UGWorkspace + UGRuname + 'timeslices.txt', "r")
    ContentUG = TimesUG.read()
    TimeslicesUG = ContentUG.split(", ")
    TimesUG.close()
    print(TimeslicesUG)
    TimesUGAD = open(UGADWorkspace + UGADRuname + 'timeslices.txt', "r")
    ContentUGAD = TimesUGAD.read()
    TimeslicesUGAD = ContentUGAD.split(", ")
    TimesUGAD.close()
    print(TimeslicesUGAD)

    TimeslicesUG[0] = 100000000
    TimeslicesUGAD[0] = 200000000

    #For each timeslice: (a)Load xyz data, (b) Calculate Mean and Max Z
    ZmeansBLOCK = []
    ZmaxsBLOCK = []
    ZmeansUG = []
    ZmaxsUG = []
    ZmeansUGAD = []
    ZmaxsUGAD = []

    for i in TimeslicesBLOCK:
        xyzBLOCK = pd.read_csv(BLOCKWorkspace + BLOCKRuname + '_time' + str(i) + '.csv')
        ymin = InitialFaultPositionBLOCK - ((int(str(i)))*ProprateBLOCK)
        ymax = ymin + WedgeWidth

        xyzUpliftZoneBLOCK = xyzBLOCK.loc[(xyzBLOCK['y'] >= ymin) & (xyzBLOCK['y'] <= ymax)]
        MeanZBLOCK = xyzUpliftZoneBLOCK['z'].mean()
        MaxZBLOCK = xyzUpliftZoneBLOCK['z'].max()
        ZmeansBLOCK.append(MeanZBLOCK)
        ZmaxsBLOCK.append(MaxZBLOCK)
        print(i)

    for i in TimeslicesUG:
        xyzUG = pd.read_csv(UGWorkspace + UGRuname + '_time' + str(i) + '.csv')
        ymin = InitialFaultPositionUG - ((int(str(i)))*ProprateUG)
        ymax = ymin + WedgeWidth

        xyzUpliftZoneUG = xyzUG.loc[(xyzUG['y'] >= ymin) & (xyzUG['y'] <= ymax)]
        MeanZUG = xyzUpliftZoneUG['z'].mean()
        MaxZUG = xyzUpliftZoneUG['z'].max()
        ZmeansUG.append(MeanZUG)
        ZmaxsUG.append(MaxZUG)
        print(i)

    for i in TimeslicesUGAD:
        xyzUGAD = pd.read_csv(UGADWorkspace + UGADRuname + '_time' + str(i) + '.csv')
        ymin = InitialFaultPositionUGAD - (((int(str(i)))-200000000)*ProprateUGAD)
        ymax = ymin + WedgeWidth

        xyzUpliftZoneUGAD = xyzUGAD.loc[(xyzUGAD['y'] >= ymin) & (xyzUGAD['y'] <= ymax)]
        MeanZUGAD = xyzUpliftZoneUGAD['z'].mean()
        MaxZUGAD = xyzUpliftZoneUGAD['z'].max()
        ZmeansUGAD.append(MeanZUGAD)
        ZmaxsUGAD.append(MaxZUGAD)
        print(i)

    #Time to Plot change in mean and maximum elevation and total material in the model by time
    TimeslicesBLOCK = [(int(i))+BLOCKTimeAdjust for i in TimeslicesBLOCK]
    TimeslicesUG = [(int(i))+UGTimeAdjust for i in TimeslicesUG]
    TimeslicesUGAD = [(int(i))+UGADTimeAdjust for i in TimeslicesUGAD]

    fig, ax1 = plt.subplots()

    ax1.plot(TimeslicesBLOCK, ZmeansBLOCK, color = 'k')
    ax1.plot(TimeslicesBLOCK, ZmaxsBLOCK, color = 'r')
    ax1.plot(TimeslicesUG, ZmeansUG, color = 'k')
    ax1.plot(TimeslicesUG, ZmaxsUG, color = 'r')
    ax1.plot(TimeslicesUGAD, ZmeansUGAD, color = 'k')
    ax1.plot(TimeslicesUGAD, ZmaxsUGAD, color = 'r')
    ax1.set_xlim(0,300000000)
    ax1.set_ylim(0,1500)
    ax1.axvspan(100000000, 200000000, color = 'lavender')
    ax1.axvspan(200000000, 300000000, color = 'lightsteelblue')
    ax2 = ax1.twinx()
    ax2.plot(TimeBLOCK, MassBLOCK, color = 'b')
    ax2.plot(TimeUG, MassUG, color = 'b')
    ax2.plot(TimeUGAD, MassUGAD, color = 'b')
    ax1.set_xlabel('Time (Myr)')
    ax1.set_ylabel('Elevation (m)')
    ax2.set_ylabel('Landscape Volume (m\u00b3)')
    ax1.set_title(FigureTitle)

    #Write the figure to file
    plt.savefig(Figure + '.png')
    print('figure saved')