
def Roundup(x):
    import math
    return int(math.ceil(x/100)) * 100

def MDD(Workspace,Runame,FaultStartingPosition,FaultPropRate,StartTime,WedgeWidth,Xmin,Xmax):
    # 1Import relevant modules, check extensions available, and create environment workspace
    import arcpy
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    arcpy.CheckOutExtension('3D')
    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("Highways")
    Workspace = arcpy.env.workspace = Workspace
    # arcpy.ImportToolbox("C:\Program Files (x86)\ArcGIS\Desktop10.6\ArcToolbox\Toolboxes\Spatial Analyst Tools.tbx")
    print('imported')

    # 2Import the data. This is the txt file of timeslices from the CHILD_to_xyz script.
    Times = open(Runame + "timeslices2run.txt", "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    # 3Calculate Slope for each timestep DEM
    ntotal = len(Timeslices)
    n = 0

    # Lists to be populated at each timestep
    MDDyMeanPerTimestep = []
    MDDyMaxPerTimestep = []
    MDDyMinPerTimestep = []

    while n < ntotal:
        # Clip DEM for 1 cell at top and then at bottom of raster (the outlet points of the raster)
        Time = Timeslices[n]
        DEM = arcpy.Raster(Workspace + str(Time) + 'DEM.tif')

        # Section 1 Identify the MDD
        # Start of MDD Identification (as origionating from MDDextraction.py script)
        DEMtopoutlet = 'DEMtopoutlet' + str(Time) + '.tif'
        DEMbottomoutlet = 'DEMbottomoutlet' + str(Time) + '.tif'

        FaultPosition = FaultStartingPosition - (FaultPropRate * ((int(Time))-StartTime))

        YminBot = Roundup(FaultPosition)
        YmaxBot = YminBot + 99
        XminBot = Xmin
        XmaxBot = Xmax
        rectangleDim = str(XminBot) + " " + str(YminBot) + " " + str(XmaxBot) + " " + str(YmaxBot)
        print(rectangleDim)
        arcpy.Clip_management(DEM, rectangle=rectangleDim, out_raster=DEMbottomoutlet)
        print(DEMbottomoutlet)

        YminTop = YminBot + (WedgeWidth - 100)
        YmaxTop = YmaxBot + (WedgeWidth - 100)
        XminTop = Xmin
        XmaxTop = Xmax
        rectangleDimTop = str(XminTop) + " " + str(YminTop) + " " + str(XmaxTop) + " " + str(YmaxTop)
        print(rectangleDimTop)
        arcpy.Clip_management(DEM, rectangle=rectangleDimTop, out_raster=DEMtopoutlet)
        print(DEMtopoutlet)

        # Convert the rasters represetning outlet nodes to points
        TopOutPoints = 'TopOutPoints' + str(Time) + '.shp'
        BotOutPoints = 'BotOutPoints' + str(Time) + '.shp'
        arcpy.RasterToPoint_conversion(DEMtopoutlet, TopOutPoints)
        print(TopOutPoints)
        arcpy.RasterToPoint_conversion(DEMbottomoutlet, BotOutPoints)
        print(BotOutPoints)

        # Calculate fill and flow direction for the DEM
        DEMFill = arcpy.sa.Fill(DEM)
        print('FillDEM')
        FD = arcpy.sa.FlowDirection(DEMFill)
        print('FD')

        # Use points to identify Upslope Contributing Areas
        TopWSs = arcpy.sa.Watershed(FD, TopOutPoints)
        print('TopWS')
        BotWSs = arcpy.sa.Watershed(FD, BotOutPoints)
        print('BotWS')

        # Convert Watershed rasters to polygons
        TopWSsPol = arcpy.RasterToPolygon_conversion(TopWSs, simplify="NO_SIMPLIFY")
        print('TopWSsPol')
        BotWSsPol = arcpy.RasterToPolygon_conversion(BotWSs, simplify="NO_SIMPLIFY")
        print('BotWSsPol')

        # Aggregate UCAs for each side of the MDD
        TopUCA = arcpy.AggregatePolygons_cartography(TopWSsPol, aggregation_distance=5)
        print('TopUCA')
        BotUCA = arcpy.AggregatePolygons_cartography(BotWSsPol, aggregation_distance=5)
        print('BotUCA')

        # Buffer Aggregated UCAs
        TopUCAbuff = arcpy.Buffer_analysis(TopUCA, buffer_distance_or_field=10)
        print('TopUCAbuff')
        BotUCAbuff = arcpy.Buffer_analysis(BotUCA, buffer_distance_or_field=10)
        print('BotUCAbuff')

        # Clip one UCAbuff using the other UCAbuff
        MDDPoly = arcpy.Clip_analysis(TopUCAbuff, BotUCAbuff, )
        print('MDDPoly')

        # MDD Polygon to Raster Conversion
        MDDRast = arcpy.PolygonToRaster_conversion(MDDPoly, "FID", cell_assignment="MAXIMUM_AREA")
        print('MDDRast')

        # MDD Raster to Polyline Conversion
        MDDPolylineName = 'MDD' + str(Time) + '.shp'
        arcpy.RasterToPolyline_conversion(MDDRast, MDDPolylineName, background_value="NODATA")
        print(MDDPolylineName)

        # MDD Raster to Point
        MDDPointName = 'MDDpoints' + str(Time) + '.shp'
        arcpy.RasterToPoint_conversion(MDDRast, MDDPointName)
        arcpy.AddXY_management(MDDPointName)
        print(MDDPointName)

        # Export Point layer as CSV file. Apply statistical analysis such that the mean y coordinate for each x coordinate value is calculated.
        GeoDateOut = 'GeoDateMDD' + str(Time) + '.gdb'
        arcpy.CreateFileGDB_management(Workspace, GeoDateOut)
        arcpy.FeatureClassToGeodatabase_conversion(MDDPointName, GeoDateOut)
        CSVfileNameMDD = 'MDD' + str(Time) + '.csv'
        InTable = str(GeoDateOut) + '/MDDpoints' + str(Time) + '_shp'
        arcpy.TableToTable_conversion(InTable, Workspace, CSVfileNameMDD)

        SACSVfileNameMDD = 'MDDsa' + str(Time) + '.csv'
        arcpy.Statistics_analysis(MDDPointName, SACSVfileNameMDD,
                                  [["POINT_Y", "MIN"], ["POINT_Y", "MAX"], ["POINT_Y", "MEAN"]], "POINT_X")
        print(SACSVfileNameMDD)

        # To prepare data for plotting MDD (distance from x axis by mean y coordinate)
        Time = Timeslices[n]
        MDDfilename = Workspace + "MDDsa" + str(Time) + ".csv"
        MDDsa = pd.read_csv(str(MDDfilename))
        # Drop X values <1000 and >4000 to limit the impact from the edge of the model
        MDDsaCent = MDDsa.drop(MDDsa.index[MDDsa['POINT_X'] < 5000])
        MDDsaCentral = MDDsaCent.drop(MDDsaCent.index[MDDsaCent['POINT_X'] > 45000])
        # Calculate Mean, Max, and Min y coordinate of the MDD and append these to a list for all timesteps
        MDDyMean = MDDsaCentral['MEAN_POINT_Y'].mean()
        MDDyMax = MDDsaCentral['MEAN_POINT_Y'].max()
        MDDyMin = MDDsaCentral['MEAN_POINT_Y'].min()
        MDDyMeanPerTimestep.append(MDDyMean)
        MDDyMaxPerTimestep.append(MDDyMax)
        MDDyMinPerTimestep.append(MDDyMin)
        print(Time)
        print('Mean' + str(MDDyMean))
        print('Max' + str(MDDyMax))
        print('Min' + str(MDDyMin))

        CSVMDDcentralName = Workspace + 'MDDcentralmean' + str(Time) + '.csv'
        MDDsaCentral.to_csv(CSVMDDcentralName)

        # Convert the central statistical analysis CSV file to point layer and then to polyline layer. This is using the mean value for each x coordinate
        XYMDDpointsName = 'MDDxyCentPoints' + str(Time) + '.shp'
        arcpy.MakeXYEventLayer_management(CSVMDDcentralName, "POINT_X", "MEAN_POINT_Y", XYMDDpointsName)
        print(XYMDDpointsName)

        MDDPoints = 'MDDpointsCent' + str(Time) + '.shp'
        arcpy.CopyFeatures_management(XYMDDpointsName, MDDPoints)

        MDDline = 'MDDCentPolyline' + str(Time) + '.shp'
        arcpy.PointsToLine_management(XYMDDpointsName, MDDline)
        print(MDDline)
        # End of MDD identification section (as origionating from MDDextraction.py script)

        n = n + 1

    # Plot MDD data from CSV files with mean y coordinates for each x coordinate, remove entries that may have been infleunced by the boundary of the model (here x<1000 or x>4000)(Mean, Min, Max)
    OutMDDyMeanFileName = Workspace + Runame + 'MDDyMeanWithTime.txt'
    MDDyMeanFileName = open(OutMDDyMeanFileName, 'w')
    MDDyMeanFileName.write(str(MDDyMeanPerTimestep))
    MDDyMeanFileName.close()

    OutMDDyMaxFileName = Workspace + Runame + 'MDDyMaxWithTime.txt'
    MDDyMaxFileName = open(OutMDDyMaxFileName, 'w')
    MDDyMaxFileName.write(str(MDDyMaxPerTimestep))
    MDDyMaxFileName.close()

    OutMDDyMinFileName = Workspace + Runame + 'MDDyMinWithTime.txt'
    MDDyMinFileName = open(OutMDDyMinFileName, 'w')
    MDDyMinFileName.write(str(MDDyMinPerTimestep))
    MDDyMinFileName.close()

def GM(Workspace, Runame, FaultStartingPosition, FaultPropRate, StartTime, WedgeWidth, Xmin, Xmax):
    # 1Import relevant modules, check extensions available, and create environment workspace
    import arcpy
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    arcpy.CheckOutExtension('3D')
    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("Highways")
    Workspace = arcpy.env.workspace = Workspace
    # arcpy.ImportToolbox("C:\Program Files (x86)\ArcGIS\Desktop10.6\ArcToolbox\Toolboxes\Spatial Analyst Tools.tbx")
    print('imported')

    # 2Import the data. This is the txt file of timeslices from the CHILD_to_xyz script.
    Times = open(Workspace + Runame + "timeslicesGM.txt", "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    # 3Calculate Slope for each timestep DEM
    ntotal = len(Timeslices)
    n = 0

    # Lists to be populated at each timestep
    MDDyMeanPerTimestepGM = []
    MDDyMaxPerTimestepGM = []
    MDDyMinPerTimestepGM = []

    while n < ntotal:
        # Clip DEM for 1 cell at top and then at bottom of raster (the outlet points of the raster)
        Time = Timeslices[n]
        DEM = arcpy.Raster(Workspace + str(Time) + 'DEM.tif')

        # Section 1 Identify the MDD
        # Start of MDD Identification (as origionating from MDDextraction.py script)
        DEMtopoutlet = 'DEMtopoutlet' + str(Time) + '.tif'
        DEMbottomoutlet = 'DEMbottomoutlet' + str(Time) + '.tif'

        FaultPosition = FaultStartingPosition - (FaultPropRate * (int(Time)-StartTime))

        YminBot = Roundup(FaultPosition)
        YmaxBot = YminBot + 99
        XminBot = Xmin
        XmaxBot = Xmax
        rectangleDim = str(XminBot) + " " + str(YminBot) + " " + str(XmaxBot) + " " + str(YmaxBot)
        print(rectangleDim)
        arcpy.Clip_management(DEM, rectangle=rectangleDim, out_raster=DEMbottomoutlet)
        print(DEMbottomoutlet)

        YminTop = YminBot + (WedgeWidth - 100)
        YmaxTop = YmaxBot + (WedgeWidth - 100)
        XminTop = Xmin
        XmaxTop = Xmax
        rectangleDimTop = str(XminTop) + " " + str(YminTop) + " " + str(XmaxTop) + " " + str(YmaxTop)
        print(rectangleDimTop)
        arcpy.Clip_management(DEM, rectangle=rectangleDimTop, out_raster=DEMtopoutlet)
        print(DEMtopoutlet)

        # Convert the rasters represetning outlet nodes to points
        TopOutPoints = 'TopOutPoints' + str(Time) + '.shp'
        BotOutPoints = 'BotOutPoints' + str(Time) + '.shp'
        arcpy.RasterToPoint_conversion(DEMtopoutlet, TopOutPoints)
        print(TopOutPoints)
        arcpy.RasterToPoint_conversion(DEMbottomoutlet, BotOutPoints)
        print(BotOutPoints)

        # Calculate fill and flow direction for the DEM
        DEMFill = arcpy.sa.Fill(DEM)
        print('FillDEM')
        FD = arcpy.sa.FlowDirection(DEMFill)
        print('FD')

        # Use points to identify Upslope Contributing Areas
        TopWSs = arcpy.sa.Watershed(FD, TopOutPoints)
        print('TopWS')
        BotWSs = arcpy.sa.Watershed(FD, BotOutPoints)
        print('BotWS')

        # Convert Watershed rasters to polygons
        TopWSsPol = arcpy.RasterToPolygon_conversion(TopWSs, simplify="NO_SIMPLIFY")
        print('TopWSsPol')
        BotWSsPol = arcpy.RasterToPolygon_conversion(BotWSs, simplify="NO_SIMPLIFY")
        print('BotWSsPol')

        # Aggregate UCAs for each side of the MDD
        TopUCA = arcpy.AggregatePolygons_cartography(TopWSsPol, aggregation_distance=5)
        print('TopUCA')
        BotUCA = arcpy.AggregatePolygons_cartography(BotWSsPol, aggregation_distance=5)
        print('BotUCA')

        # Buffer Aggregated UCAs
        TopUCAbuff = arcpy.Buffer_analysis(TopUCA, buffer_distance_or_field=10)
        print('TopUCAbuff')
        BotUCAbuff = arcpy.Buffer_analysis(BotUCA, buffer_distance_or_field=10)
        print('BotUCAbuff')

        # Clip one UCAbuff using the other UCAbuff
        MDDPoly = arcpy.Clip_analysis(TopUCAbuff, BotUCAbuff, )
        print('MDDPoly')

        # MDD Polygon to Raster Conversion
        MDDRast = arcpy.PolygonToRaster_conversion(MDDPoly, "FID", cell_assignment="MAXIMUM_AREA")
        print('MDDRast')

        # MDD Raster to Polyline Conversion
        MDDPolylineName = 'MDD' + str(Time) + '.shp'
        arcpy.RasterToPolyline_conversion(MDDRast, MDDPolylineName, background_value="NODATA")
        print(MDDPolylineName)

        # MDD Raster to Point
        MDDPointName = 'MDDpoints' + str(Time) + '.shp'
        arcpy.RasterToPoint_conversion(MDDRast, MDDPointName)
        arcpy.AddXY_management(MDDPointName)
        print(MDDPointName)

        # Export Point layer as CSV file. Apply statistical analysis such that the mean y coordinate for each x coordinate value is calculated.
        GeoDateOut = 'GeoDateMDD' + str(Time) + '.gdb'
        arcpy.CreateFileGDB_management(Workspace, GeoDateOut)
        arcpy.FeatureClassToGeodatabase_conversion(MDDPointName, GeoDateOut)
        CSVfileNameMDD = 'MDD' + str(Time) + '.csv'
        InTable = 'MDDpoints' + str(Time) + '.shp'
        arcpy.TableToTable_conversion(InTable, Workspace, CSVfileNameMDD)

        SACSVfileNameMDD = 'MDDsa' + str(Time) + '.csv'
        arcpy.Statistics_analysis(MDDPointName, SACSVfileNameMDD, [["POINT_Y", "MIN"], ["POINT_Y", "MAX"], ["POINT_Y", "MEAN"]], "POINT_X")
        print(SACSVfileNameMDD)

        # To prepare data for plotting MDD (distance from x axis by mean y coordinate)
        Time = Timeslices[n]
        MDDfilename = Workspace + "MDDsa" + str(Time) + ".csv"
        MDDsa = pd.read_csv(str(MDDfilename))
        # Drop X values <1000 and >4000 to limit the impact from the edge of the model
        MDDsaCent = MDDsa.drop(MDDsa.index[MDDsa['POINT_X'] < 5000])
        MDDsaCentral = MDDsaCent.drop(MDDsaCent.index[MDDsaCent['POINT_X'] > 45000])
        # Calculate Mean, Max, and Min y coordinate of the MDD and append these to a list for all timesteps
        MDDyMean = MDDsaCentral['MEAN_POINT_Y'].mean()
        MDDyMax = MDDsaCentral['MEAN_POINT_Y'].max()
        MDDyMin = MDDsaCentral['MEAN_POINT_Y'].min()
        MDDyMeanPerTimestepGM.append(MDDyMean)
        MDDyMaxPerTimestepGM.append(MDDyMax)
        MDDyMinPerTimestepGM.append(MDDyMin)
        print(Time)

        CSVMDDcentralName = Workspace + 'MDDcentralmean' + str(Time) + '.csv'
        MDDsaCentral.to_csv(CSVMDDcentralName)

        # Convert the central statistical analysis CSV file to point layer and then to polyline layer. This is using the mean value for each x coordinate
        XYMDDpointsName = 'MDDxyCentPoints' + str(Time) + '.shp'
        arcpy.MakeXYEventLayer_management(CSVMDDcentralName, "POINT_X", "MEAN_POINT_Y", XYMDDpointsName)
        print(XYMDDpointsName)

        MDDPoints = 'MDDpointsCent' + str(Time) + '.shp'
        arcpy.CopyFeatures_management(XYMDDpointsName, MDDPoints)

        MDDline = 'MDDCentPolyline' + str(Time) + '.shp'
        arcpy.PointsToLine_management(XYMDDpointsName, MDDline)
        print(MDDline)
        # End of MDD identification section (as origionating from MDDextraction.py script)

        FA = arcpy.sa.FlowAccumulation(FD)
        print('FA')
        UCA100000 = FA >= 10  # 40 if using DEM of resolution 50 m. 10 if using DEM of resolution 100 m.
        FSUCA = arcpy.sa.FocalStatistics(UCA100000, arcpy.sa.NbrRectangle(3, 3, "CELL"), "SUM")
        Conhead = arcpy.sa.Con(FSUCA, UCA100000, "", "Value=2")
        HeadOnly = arcpy.sa.ExtractByAttributes(Conhead, "Value=1")
        HeadPoints = 'HeadPoints' + str(Time) + '.shp'
        arcpy.RasterToPoint_conversion(HeadOnly, HeadPoints, "Value")
        print(HeadPoints)

        # Identify watersheds from headpoints, polygonise watersheds find centre points and only take those less than 1000 m from the MDD and intersecting the MDD
        Watersheds = arcpy.sa.Watershed(FD, HeadPoints, "pointid")
        WatershedsPoly = 'WatershedsPoly' + str(Time) + '.shp'
        arcpy.RasterToPolygon_conversion(Watersheds, WatershedsPoly, "NO_SIMPLIFY")
        print(WatershedsPoly)
        WatershedsPoint = 'WatershedsPoint' + str(Time) + '.shp'
        arcpy.FeatureToPoint_management(WatershedsPoly, WatershedsPoint, "CENTROID")
        print(WatershedsPoint)

        arcpy.Near_analysis(WatershedsPoint, MDDline)
        WatershedsPointFL = 'WatershedsPointFL' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(WatershedsPoint, WatershedsPointFL)
        arcpy.SelectLayerByAttribute_management(WatershedsPointFL, "", '"NEAR_DIST"<1000')
        WatershedPointsNearMDD = 'WatershedPointsNearMDD' + str(Time) + '.shp'
        arcpy.CopyFeatures_management(WatershedsPointFL, WatershedPointsNearMDD)

        WatershedsPolyFL = 'WatershedsPolyFL' + str(Time) + '.shp'
        WatershedPointsNearMDDFL = 'WatershedPointsNearMDDFL' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(WatershedsPoly, WatershedsPolyFL)
        arcpy.MakeFeatureLayer_management(WatershedPointsNearMDD, WatershedPointsNearMDDFL)
        arcpy.SelectLayerByLocation_management(WatershedsPolyFL, "", WatershedPointsNearMDDFL)
        WatershedsPolyNearMDD = 'WatershedsPolyNearMDD' + str(Time) + '.shp'
        arcpy.CopyFeatures_management(WatershedsPolyFL, WatershedsPolyNearMDD)
        print(WatershedsPolyNearMDD)

        # Near distance between polygons and the MDD to select only headwaters that border (within 2 cells) the MDD
        arcpy.Near_analysis(WatershedsPolyNearMDD, MDDline)
        WatershedsPolyNearMDDFL = 'WatershedsPolyNearMDDFL' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(WatershedsPolyNearMDD, WatershedsPolyNearMDDFL)
        arcpy.SelectLayerByAttribute_management(WatershedsPolyNearMDDFL, "", '"NEAR_DIST"<100')
        Headwaters = 'Headwaters' + str(Time) + '.shp'
        arcpy.CopyFeatures_management(WatershedsPolyNearMDDFL, Headwaters)
        print(Headwaters)

        # Calculate Slope and Relief Rasters. No need to clip these as across divide analysis only applied to the central parts of the rasters
        # For Slope
        Slope = 'Slope' + str(Time) + '.tif'
        arcpy.Slope_3d(DEM, Slope)
        # For Relief
        # For Relief
        nbr = arcpy.sa.NbrCircle(4, "CELL")
        FSmin = arcpy.sa.FocalStatistics(DEM, nbr, "MINIMUM", "DATA")
        FSmax = arcpy.sa.FocalStatistics(DEM, nbr, "MAXIMUM", "DATA")
        Relief4Cell = FSmax - FSmin
        ReliefOutName = 'Relief4c' + str(Time) + '.tif'
        Relief4Cell.save(ReliefOutName)
        print(ReliefOutName)

        # Calculate Zonal Statistics for the headwaters bordering the MDD (Mean Z, S, R)
        HeadwatersZ = 'HeadwatersZ' + str(Time) + '.shp'
        HeadwatersS = 'HeadwatersS' + str(Time) + '.shp'
        HeadwatersR = 'HeadwatersR' + str(Time) + '.shp'
        arcpy.CopyFeatures_management(Headwaters, HeadwatersZ)
        arcpy.CopyFeatures_management(Headwaters, HeadwatersS)
        arcpy.CopyFeatures_management(Headwaters, HeadwatersR)
        ZSZ = arcpy.sa.ZonalStatistics(HeadwatersZ, "FID", DEMFill, "MEAN")
        print(HeadwatersZ)
        ZSS = arcpy.sa.ZonalStatistics(HeadwatersS, "FID", Slope, "MEAN")
        print(HeadwatersS)
        ZSR = arcpy.sa.ZonalStatistics(HeadwatersR, "FID", Relief4Cell, "MEAN")
        print(HeadwatersR)

        # Extract the Zonal Statistic means to point for the headwaters
        HeadwaterPoints = 'HeadwaterPoints' + str(Time) + '.shp'
        arcpy.FeatureToPoint_management(Headwaters, HeadwaterPoints)
        HeadwaterPointsZ = 'HeadwaterPointsZ' + str(Time) + '.shp'
        HeadwaterPointsS = 'HeadwaterPointsS' + str(Time) + '.shp'
        HeadwaterPointsR = 'HeadwaterPointsR' + str(Time) + '.shp'
        arcpy.sa.ExtractValuesToPoints(HeadwaterPoints, ZSZ, HeadwaterPointsZ)
        print(HeadwaterPointsZ)
        arcpy.sa.ExtractValuesToPoints(HeadwaterPoints, ZSS, HeadwaterPointsS)
        print(HeadwaterPointsS)
        arcpy.sa.ExtractValuesToPoints(HeadwaterPoints, ZSR, HeadwaterPointsR)
        print(HeadwaterPointsR)

        # Split Headwater points based on side of the MDD
        # Layers to be created
        HeadwaterPointsZtop = 'HeadwaterPointsZtop' + str(Time) + '.shp'
        HeadwaterPointsZtop2 = 'HeadwaterPointsZtop2' + str(Time) + '.shp'
        HeadwaterPointsZbot = 'HeadwaterPointsZbot' + str(Time) + '.shp'
        HeadwaterPointsStop = 'HeadwaterPointsStop' + str(Time) + '.shp'
        HeadwaterPointsStop2 = 'HeadwaterPointsStop2' + str(Time) + '.shp'
        HeadwaterPointsSbot = 'HeadwaterPointsSbot' + str(Time) + '.shp'
        HeadwaterPointsRtop = 'HeadwaterPointsRtop' + str(Time) + '.shp'
        HeadwaterPointsRtop2 = 'HeadwaterPointsRtop2' + str(Time) + '.shp'
        HeadwaterPointsRbot = 'HeadwaterPointsRbot' + str(Time) + '.shp'

        # Set up Feature Layers
        TopUCAFL = 'TopUCAFL' + str(Time) + '.shp'
        BotUCAFL = 'BotUCAFL' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(TopUCA, TopUCAFL)
        arcpy.MakeFeatureLayer_management(BotUCA, BotUCAFL)
        HeadwaterPointsZFL = 'HeadwaterPointsZFL' + str(Time) + '.shp'
        HeadwaterPointsSFL = 'HeadwaterPointsSFL' + str(Time) + '.shp'
        HeadwaterPointsRFL = 'HeadwaterPointsRFL' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsZ, HeadwaterPointsZFL)
        arcpy.MakeFeatureLayer_management(HeadwaterPointsS, HeadwaterPointsSFL)
        arcpy.MakeFeatureLayer_management(HeadwaterPointsR, HeadwaterPointsRFL)
        print('Feature Layers set up for Gilbert Metric Analysis')

        # Select headwater points on each side of the MDD for Z,S,R and export to new layer
        arcpy.SelectLayerByLocation_management(HeadwaterPointsZFL, "intersect", TopUCAFL)
        arcpy.CopyFeatures_management(HeadwaterPointsZFL, HeadwaterPointsZtop)
        arcpy.CopyFeatures_management(HeadwaterPointsZFL, HeadwaterPointsZtop2)
        arcpy.SelectLayerByLocation_management(HeadwaterPointsZFL, "intersect", BotUCAFL)
        arcpy.CopyFeatures_management(HeadwaterPointsZFL, HeadwaterPointsZbot)

        arcpy.SelectLayerByLocation_management(HeadwaterPointsSFL, "intersect", TopUCAFL)
        arcpy.CopyFeatures_management(HeadwaterPointsSFL, HeadwaterPointsStop)
        arcpy.CopyFeatures_management(HeadwaterPointsSFL, HeadwaterPointsStop2)
        arcpy.SelectLayerByLocation_management(HeadwaterPointsSFL, "intersect", BotUCAFL)
        arcpy.CopyFeatures_management(HeadwaterPointsSFL, HeadwaterPointsSbot)

        arcpy.SelectLayerByLocation_management(HeadwaterPointsRFL, "intersect", TopUCAFL)
        arcpy.CopyFeatures_management(HeadwaterPointsRFL, HeadwaterPointsRtop)
        arcpy.CopyFeatures_management(HeadwaterPointsRFL, HeadwaterPointsRtop2)
        arcpy.SelectLayerByLocation_management(HeadwaterPointsRFL, "intersect", BotUCAFL)
        arcpy.CopyFeatures_management(HeadwaterPointsRFL, HeadwaterPointsRbot)
        print('Headwater Points with Z,S,R prepared for Gilbert Metric analysis')

        # Near
        arcpy.Near_analysis(HeadwaterPointsZtop2, HeadwaterPointsZbot)
        arcpy.Near_analysis(HeadwaterPointsZbot, HeadwaterPointsZtop)
        arcpy.Near_analysis(HeadwaterPointsStop2, HeadwaterPointsSbot)
        arcpy.Near_analysis(HeadwaterPointsSbot, HeadwaterPointsStop)
        arcpy.Near_analysis(HeadwaterPointsRtop2, HeadwaterPointsRbot)
        arcpy.Near_analysis(HeadwaterPointsRbot, HeadwaterPointsRtop)
        print('Near distances for neighbouring across divide headwaters calculated')

        # Make Feature Layers for Joining Tables
        ZbasinsTop2 = 'ZbasinsTop2' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsZtop2, ZbasinsTop2)
        ZbasinsBot = 'ZbasinsBot' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsZbot, ZbasinsBot)
        ZbasinsTop = 'ZbasinsTop' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsZtop, ZbasinsTop)
        SbasinsTop2 = 'SbasinsTop2' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsStop2, SbasinsTop2)
        SbasinsBot = 'SbasinsBot' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsSbot, SbasinsBot)
        SbasinsTop = 'SbasinsTop' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsStop, SbasinsTop)
        RbasinsTop2 = 'RbasinsTop2' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsRtop2, RbasinsTop2)
        RbasinsBot = 'RbasinsBot' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsRbot, RbasinsBot)
        RbasinsTop = 'RbasinsTop' + str(Time) + '.shp'
        arcpy.MakeFeatureLayer_management(HeadwaterPointsRtop, RbasinsTop)

        # Join Tables
        arcpy.AddJoin_management(ZbasinsTop2, "NEAR_FID", ZbasinsBot, "FID")
        arcpy.AddJoin_management(ZbasinsBot, "NEAR_FID", ZbasinsTop, "FID")
        arcpy.AddJoin_management(SbasinsTop2, "NEAR_FID", SbasinsBot, "FID")
        arcpy.AddJoin_management(SbasinsBot, "NEAR_FID", SbasinsTop, "FID")
        arcpy.AddJoin_management(RbasinsTop2, "NEAR_FID", RbasinsBot, "FID")
        arcpy.AddJoin_management(RbasinsBot, "NEAR_FID", RbasinsTop, "FID")
        print('Headwater Pairs Identified for all near MDD headwaters')

        # Export Tables as Csv files
        # TopBot files; first metric entry relates to the basin to the top of the MMD, the second relates to the basin to the bottom of the MDD
        # BotTop files; first metric entry relates to the basin to the bottom of the MDD, the seconf relates to the basin to the top of the MDD
        ZjoinedTopBot = 'ZjoinedTopBot' + str(Time) + '.csv'
        arcpy.TableToTable_conversion(ZbasinsTop2, Workspace, ZjoinedTopBot)
        print(ZjoinedTopBot)
        ZjoinedBotTop = 'ZjoinedBotTop' + str(Time) + '.csv'
        arcpy.TableToTable_conversion(ZbasinsBot, Workspace, ZjoinedBotTop)
        print(ZjoinedBotTop)
        SjoinedTopBot = 'SjoinedTopBot' + str(Time) + '.csv'
        arcpy.TableToTable_conversion(SbasinsTop2, Workspace, SjoinedTopBot)
        print(SjoinedTopBot)
        SjoinedBotTop = 'SjoinedBotTop' + str(Time) + '.csv'
        arcpy.TableToTable_conversion(SbasinsBot, Workspace, SjoinedBotTop)
        print(SjoinedBotTop)
        RjoinedTopBot = 'RjoinedTopBot' + str(Time) + '.csv'
        arcpy.TableToTable_conversion(RbasinsTop2, Workspace, RjoinedTopBot)
        print(RjoinedTopBot)
        RjoinedBotTop = 'RjoinedBotTop' + str(Time) + '.csv'
        arcpy.TableToTable_conversion(RbasinsBot, Workspace, RjoinedBotTop)
        print(RjoinedBotTop)
        print('Headwater Pair Tables exported to csv files')

        n = n + 1

def MDDplot(MDDlocationTable, Figure, FigureTitle):
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    MDDLocationTable = pd.read_csv(MDDlocationTable)

    fig1, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1, 3]})

    Time = MDDLocationTable["Time"] / 1000000
    MDDmean = MDDLocationTable["MeanMDDrelFault"] / 1000
    MDDmin = MDDLocationTable["MinMDDrelFault"] / 1000
    MDDmax = MDDLocationTable["MaxMDDrelFault"] / 1000
    MDDspan = MDDLocationTable["MDDspan"] / 1000

    ax2.plot(Time, MDDmean, color='black')
    ax2.plot(Time, MDDmin, color='blue')
    ax2.plot(Time, MDDmax, color='red')
    ax1.plot(Time, MDDspan, color='black')
    ax1.axvspan(100, 200, color='lavender')
    ax1.axvspan(200, 300, color='lightsteelblue')
    ax2.axvspan(100, 200, color='lavender')
    ax2.axvspan(200, 300, color='lightsteelblue')
    ax2.set_xlim(0, 300)
    ax1.set_xlim(0, 300)
    ax1.set_ylim(0, 10)
    ax2.set_ylim(0, 20)
    ax2.set_ylabel('MDD Location Relative To Fault (km)')
    ax1.set_ylabel('MDD Span (km)')
    ax2.set_xlabel('Time (Myr)')
    ax1.set_title(FigureTitle)

    #Write the figure to file
    plt.savefig(Figure + '.png')
    print('figure saved')

def GMplot(Workspace, Runame):
    import numpy as np
    from matplotlib import pyplot as plt
    import pandas as pd

    Times = open(Workspace + Runame + "timeslicesGM.txt", "r")
    content = Times.read()
    Timeslices = content.split(", ")  # [1:-1]
    Times.close()

    print(Timeslices)

    n = 0
    ntotal = len(Timeslices)

    while n < ntotal:
        Time = Timeslices[n]

        ZjoinedTopBot = 'ZjoinedTopBot' + str(Time) + '.csv'
        ZjoinedBotTop = 'ZjoinedBotTop' + str(Time) + '.csv'
        SjoinedTopBot = 'SjoinedTopBot' + str(Time) + '.csv'
        SjoinedBotTop = 'SjoinedBotTop' + str(Time) + '.csv'
        RjoinedTopBot = 'RjoinedTopBot' + str(Time) + '.csv'
        RjoinedBotTop = 'RjoinedBotTop' + str(Time) + '.csv'

        # Calculate Delta for each Gilbert Metric
        # Here it is hypothesed the divide migrates to the top of the model so:
        # Z is calculated as Top - Bot. Positive values therefore support the hypothesis
        # S is calculated as Bot - Top. Positive values therefore support the hypothesis
        # R is calculated as Bot - Top. Positive values therefore support the hypothesis
        Ztopbot = pd.read_csv(Workspace + ZjoinedTopBot)
        Zbottop = pd.read_csv(Workspace + ZjoinedBotTop)
        Stopbot = pd.read_csv(Workspace + SjoinedTopBot)
        Sbottop = pd.read_csv(Workspace + SjoinedBotTop)
        Rtopbot = pd.read_csv(Workspace + RjoinedTopBot)
        Rbottop = pd.read_csv(Workspace + RjoinedBotTop)

        HeadwaterPointsZbot0_RASTERVALU = 'HeadwaterPointsZbot' + str(Time) + '_RASTERVALU'
        HeadwaterPointsZtop0_RASTERVALU = 'HeadwaterPointsZtop' + str(Time) + '_RASTERVALU'
        HeadwaterPointsSbot0_RASTERVALU = 'HeadwaterPointsSbot' + str(Time) + '_RASTERVALU'
        HeadwaterPointsStop0_RASTERVALU = 'HeadwaterPointsStop' + str(Time) + '_RASTERVALU'
        HeadwaterPointsRbot0_RASTERVALU = 'HeadwaterPointsRbot' + str(Time) + '_RASTERVALU'
        HeadwaterPointsRtop0_RASTERVALU = 'HeadwaterPointsRtop' + str(Time) + '_RASTERVALU'

        Ztopbot['DeltaZ'] = Ztopbot['RASTERVALU'] - Ztopbot[HeadwaterPointsZbot0_RASTERVALU]
        Zbottop['DeltaZ'] = Zbottop[HeadwaterPointsZtop0_RASTERVALU] - Zbottop['RASTERVALU']
        ZtopbotDeltaList = Ztopbot['DeltaZ'].tolist()
        # print(ZtopbotDeltaList)
        ZbottopDeltaList = Zbottop['DeltaZ'].tolist()
        # print(ZbottopDeltaList)
        DeltaZ = []
        DeltaZ.extend(ZtopbotDeltaList)
        DeltaZ.extend(ZbottopDeltaList)
        print(DeltaZ)

        Stopbot['DeltaS'] = Stopbot[HeadwaterPointsSbot0_RASTERVALU] - Stopbot['RASTERVALU']
        Sbottop['DeltaS'] = Sbottop['RASTERVALU'] - Sbottop[HeadwaterPointsStop0_RASTERVALU]
        StopbotDeltaList = Stopbot['DeltaS'].tolist()
        # print(StopbotDeltaList)
        SbottopDeltaList = Sbottop['DeltaS'].tolist()
        # print(SbottopDeltaList)
        DeltaS = []
        DeltaS.extend(StopbotDeltaList)
        DeltaS.extend(SbottopDeltaList)
        print(DeltaS)

        Rtopbot['DeltaR'] = Rtopbot[HeadwaterPointsRbot0_RASTERVALU] - Rtopbot['RASTERVALU']
        Rbottop['DeltaR'] = Rbottop['RASTERVALU'] - Rbottop[HeadwaterPointsRtop0_RASTERVALU]
        RtopbotDeltaList = Rtopbot['DeltaR'].tolist()
        # print(RtopbotDeltaList)
        RbottopDeltaList = Rbottop['DeltaR'].tolist()
        # print(RbottopDeltaList)
        DeltaR = []
        DeltaR.extend(RtopbotDeltaList)
        DeltaR.extend(RbottopDeltaList)
        print(DeltaR)
        print('Gilbert Metrics calculated for all headwater pairs' + str(Time))

        # Remove duplicate pairs to only show unique pairs
        DeltaZ = list(dict.fromkeys(DeltaZ))
        DeltaS = list(dict.fromkeys(DeltaS))
        DeltaR = list(dict.fromkeys(DeltaR))

        #Save data to file
        ZfileName = Workspace + Time + 'DeltaZvalues.txt'
        SfileName = Workspace + Time + 'DeltaSvalues.txt'
        RfileName = Workspace + Time + 'DeltaRvalues.txt'

        open_fileZ = open(ZfileName, 'w')
        open_fileZ.write(", ".join( repr(e) for e in DeltaZ))
        open_fileZ.close()

        open_fileS = open(SfileName, 'w')
        open_fileS.write(", ".join( repr(e) for e in DeltaS))
        open_fileS.close()

        open_fileR = open(RfileName, 'w')
        open_fileR.write(", ".join( repr(e) for e in DeltaR))
        open_fileR.close()

        # Plot Histogram and Boxplots for Gilbert Metrics and Export to file
        No_change = 0

        Zbins = [-400, -350, -300, -250, -200, -150, -100, -50, 0, 50, 100, 150, 200, 250, 300, 350, 400]
        Sbins = [-18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18]
        Rbins = [-140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140]

        fig1, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1, 3]})

        ax2.hist(DeltaZ, bins=Zbins, color='grey', edgecolor='black', zorder=2)
        ax2.set_xlabel('\u0394 Z (m)')
        ax2.set_ylabel('Number of headwater comparisons')
        ax2.axvline(No_change, color='black')
        ax2.axis([-400, 400, 0, 80])
        ax1.boxplot(x=DeltaZ, vert=False, whis=[5, 95])
        ax1.axis([-400, 400, 0.5, 1.5])
        ax1.axvline(No_change, color='black')
        ax1.axvspan(-400, 0, color='darkblue', alpha=0.2, lw=0, zorder=1)
        ax1.axvspan(0, 400, color='gold', alpha=0.2, lw=0, zorder=1)
        ax2.axvspan(-400, 0, color='darkblue', alpha=0.2, lw=0, zorder=1)
        ax2.axvspan(0, 400, color='gold', alpha=0.2, lw=0, zorder=1)
        fig1.suptitle('Time = ' + Time)

        # Save plot to file as .png.
        plotname = Runame + 'Zdelta' + str(Time) + '.png'
        plt.savefig(plotname, dpi=72)

        fig2, (ax3, ax4) = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1, 3]})

        ax4.hist(DeltaS, bins=Sbins, color='grey', edgecolor='black', zorder=2)
        ax4.set_xlabel('\u0394 S (\u00b0)')
        ax4.set_ylabel('Number of headwater comparisons')
        ax4.axvline(No_change, color='black')
        ax4.axis([-20, 20, 0, 40])
        ax3.boxplot(x=DeltaS, vert=False, whis=[5, 95])
        ax3.axvline(No_change, color='black')
        ax3.axvspan(-400, 0, color='darkblue', alpha=0.2, lw=0, zorder=1)
        ax3.axvspan(0, 400, color='gold', alpha=0.2, lw=0, zorder=1)
        ax4.axvspan(-400, 0, color='darkblue', alpha=0.2, lw=0, zorder=1)
        ax4.axvspan(0, 400, color='gold', alpha=0.2, lw=0, zorder=1)
        fig2.suptitle('Time = ' + Time)

        # Save plot to file as .png.
        plotname = Runame + 'Sdelta' + str(Time) + '.png'
        plt.savefig(plotname, dpi=72)

        fig3, (ax5, ax6) = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios': [1, 3]})

        ax6.hist(DeltaR, bins=Rbins, color='grey', edgecolor='black', zorder=2)
        ax6.set_xlabel('\u0394 R (m)')
        ax6.set_ylabel('Number of headwater comparisons')
        ax6.axvline(No_change, color='black')
        ax6.axis([-150, 150, 0, 50])
        ax5.boxplot(x=DeltaR, vert=False, whis=[5, 95])
        ax5.axvline(No_change, color='black')
        ax5.axvspan(-400, 0, color='darkblue', alpha=0.2, lw=0, zorder=1)
        ax5.axvspan(0, 400, color='gold', alpha=0.2, lw=0, zorder=1)
        ax6.axvspan(-400, 0, color='darkblue', alpha=0.2, lw=0, zorder=1)
        ax6.axvspan(0, 400, color='gold', alpha=0.2, lw=0, zorder=1)
        fig3.suptitle('Time = ' + Time)

        # Save plot to file as .png.
        plotname = Runame + 'Rdelta' + str(Time) + '.png'
        plt.savefig(Workspace + plotname, dpi=72)

        n = n + 1

#def GMplotTimeseries(Workspace, Runame):
    #This module only exists in manual form at the moment see python file 'GMtimseriesboxplots' in 128 GM on R drive



