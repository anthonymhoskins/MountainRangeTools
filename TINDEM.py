

def TINdem(Workspace,Runame):

    # 1Import Relevant Modules, check extensions available, and create environment workspace
    import arcpy

    arcpy.CheckOutExtension('3D')
    arcpy.env.workspace = Workspace

    # 2Import the data. This is the csv xyz file (from excel conversion after CHILD_to_xyz)
    Times = open(Runame + 'timeslices.txt', "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    # 3 For each timesclice; (a)Load xy data, (b)create a TIN, and (c)create a DEM
    for i in Timeslices:
        csvData = Runame + '_time' + str(i) + '.csv'
        # (a) Load xy data (MakeXYEventLayer equivilant of Display XY Data)
        XYPointsName = str(i) + 'points'
        xypoints = arcpy.MakeXYEventLayer_management(csvData, "x", "y", XYPointsName)
        # (b) Create a TIN
        TINName = str(i) + "TIN"
        arcpy.CreateTin_3d(TINName, in_features=xypoints)
        print(TINName)
        # (c) Create a DEM
        DEMRasterName = str(i) + "DEM" + ".tif"
        arcpy.TinRaster_3d(TINName, DEMRasterName, sample_distance="CELLSIZE 100")
        print(DEMRasterName)
