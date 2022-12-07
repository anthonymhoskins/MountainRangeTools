def ERcatAv(Workspace, DEMfileName, ERxyzfileName, MDDcentralmeanfileName, LSDcatchmentsfileName):
    #Workspace is the location where files are retrieved and saved
    #DEMfileName is the file name of the DEM that this analysis is going to be applied to
    #ERxyzfileName is the file name of the ERxyz, output from ER.py
    #MDDcentralmeanfileName is the file name of the mean MDD locations associated with this timestep, output from MDDandGMextraction.py
    #LSDcatchmentsfileName is a .bil file of the catchments as extracted in LSDtopotools, using scripts 1ChannelExtraction.driver and 2CatchmentExtraction.driver

    # 1Import relevant modules, check extensions available, and create environment workspace
    import arcpy
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    arcpy.CheckOutExtension('3D')
    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("Highways")
    Workspace = arcpy.env.workspace = Workspace #i.e, 'E:/ERcatavWork/AD01/'
    # arcpy.ImportToolbox("C:\Program Files (x86)\ArcGIS\Desktop10.6\ArcToolbox\Toolboxes\Spatial Analyst Tools.tbx")
    print('imported')

    # 2 Import the data to be used; A DEM, ERxyz data (with voroni area added), and MDDcentralmean file
    DEMfile = Workspace + DEMfileName #i.e, "E:/ERcatavWork/AD01/295000000DEM.tif"
    ERxyzfile = Workspace + ERxyzfileName  #i.e, "E:/ERcatavWork/AD01/ERxyz295000000.csv"
    MDDcentralmeanfile = Workspace + MDDcentralmeanfileName  #i.e, "E:/ERcatavWork/AD01/MDDcentralmean295000000.csv"
    LSDcatchments = Workspace + LSDcatchmentsfileName #i.e, "E:/ERcatavWork/AD01/295Myrp_AllBasins.bil"

    # 3 Take the ERxyzArea csv file and make a shapefile point layer from these
    ERxyzEventLayerName = 'ERxyzEL'
    ERxyzAreaPointsSaveLoc = Workspace + 'ERxyzAreaPoints' + '.lyr'
    arcpy.MakeXYEventLayer_management(ERxyzfile, "x", "y", ERxyzEventLayerName)
    arcpy.SaveToLayerFile_management(ERxyzEventLayerName, ERxyzAreaPointsSaveLoc)
    ERxyzAPointShpName = "ERxyzAPoints.shp"
    arcpy.CopyFeatures_management(ERxyzAreaPointsSaveLoc, ERxyzAPointShpName)
    print('ERxyzAPoints.shp exported')

    # Using the LSD catchment Raster output we will clip the filled DEM to the size of the catchments
    DEMfill = arcpy.sa.Fill(DEMfile)
    print('DEM filled')
    CatPoly = Workspace + 'LSDcats.shp'
    LSDcatchments = arcpy.sa.Int(LSDcatchments)
    arcpy.RasterToPolygon_conversion(LSDcatchments, CatPoly)
    print('LSD catchments exported as shp')
    DEMlsd = Workspace + 'DEMlsd.tif'
    arcpy.Clip_management(DEMfill, rectangle=None, out_raster=DEMlsd, in_template_dataset=CatPoly,
                          clipping_geometry="ClippingGeometry")
    # arcpy.Clip_management(DEMfill,CatPoly,DEMlsd,clipping_geometry="ClippingGeometry")
    print('DEM clipped to LSD catchments')

    # 4 Using a DEM that has been filled and clipped to the catchments identified through LSDTopoTools, the following layers are created;
    # 1.Flow direction raster, 2. Flow accumulation raster, 3. Flow length downstream Raster,
    # 4. Flow accumulation weighted by sum of ER in ERxyzAPoints.shp, 5. Flow accumulation weighted by sum of Area in ERxyzAPoints.shp

    FD = arcpy.sa.FlowDirection(DEMlsd)
    print('FD created')

    FA = arcpy.sa.FlowAccumulation(FD)
    print('FA created')

    FL = arcpy.sa.FlowLength(FD, "DOWNSTREAM")
    print('FL created')

    # Rasterise ERxyzAPoints with ER as value
    ERrastOut = Workspace + "ERrast.tif"
    arcpy.env.snapRaster = DEMlsd
    arcpy.PointToRaster_conversion(ERxyzAPointShpName, "ER", ERrastOut, "MEAN", cellsize=DEMlsd)
    print('ER raster created')

    # Rasterise ERxyzAPoints with Area as value
    ArearastOut = Workspace + "Arearast.tif"
    arcpy.env.snapRaster = DEMlsd
    arcpy.PointToRaster_conversion(ERxyzAPointShpName, "Area", ArearastOut, "SUM", cellsize=DEMlsd)
    print('Area raster created')

    # To get the volume of material eroded multiply ER raster by Area Raster
    Raster1 = arcpy.Raster(ERrastOut)
    Raster2 = arcpy.Raster(ArearastOut)
    Volume = Raster1 * Raster2
    VolRastName = 'VolumeRast.tif'
    Volume.save(Workspace + VolRastName)
    print('Volume raster created')

    FAforVol = arcpy.sa.FlowAccumulation(FD, VolRastName)
    FAforVol.save(Workspace + 'VolumeFA.tif')
    print('Flow Accumulation Volume raster created')

    FAmm = 'FAm2.tif'
    UCA = ((FA) * (100 * 100))
    UCA.save(Workspace + FAmm)

    # Now we can move on to extract river channels and export raster cells representing river channels to a point shp file
    RiverRast = arcpy.sa.Con(FA, 1, where_clause="VALUE > 10")
    print('River raster created')
    RiverPoint = Workspace + 'RiverPoints.shp'
    arcpy.RasterToPoint_conversion(RiverRast, RiverPoint)
    print('River Points created')

    # Now we will extract the values from each of the rasters we are interested in to the river points shapefile

    arcpy.sa.ExtractMultiValuesToPoints(RiverPoint, [[FA, "FA"], [DEMlsd, "Z"], [FL, "FL"], [FAforVol, "ERvol"],
                                                     [LSDcatchments, "CatNum"], [FAmm, "FAm2"]])
    print('Values extracted to river points')

    # Now we can look at the MDD and calculate the distance between the MDD and the river points
    MDDEventLayerName = 'MDDEL'
    MDDSaveLoc = Workspace + 'MDDPoints' + '.lyr'
    arcpy.MakeXYEventLayer_management(MDDcentralmeanfile, "POINT_X", "MEAN_POINT_Y", MDDEventLayerName)
    arcpy.SaveToLayerFile_management(MDDEventLayerName, MDDSaveLoc)
    MDDPointShpName = Workspace + "MDDPoints.shp"
    arcpy.CopyFeatures_management(MDDSaveLoc, MDDPointShpName)
    print('MDDPoints.shp exported')

    arcpy.Near_analysis(RiverPoint, MDDPointShpName, angle="ANGLE")
    print('Distance to MDD extracted to River points')

def ERcatAvPlot(Workspace, RiverPointsFileName, ProxCatNo, DistCatNo):
    #Workspace is the workspace, i.e. "E:/ERcatavWork/
    #RiverPointsFileName is the name of the river points .txt file holding the ER data (output from ERcatAv).
    #ProxCatNo the catchment number associated with the proximal flank catchment we wish to plot data for
    #DistCatNo the catchment number associated with the distal flank catchment we wish to polot data for

    #This script will plot 4 graphs;
    #(1) Flow dist to outlet by Z with c as ERcatAv
    #(2) Distance from the MDD by Z with c as ERcatAv
    #(3) Flow dist to outlet by ERcatAv with c as ERcatAv
    #(4) Distance from the MDD by ERcatAv with c as ERcatAv
    #These will be saved to the workspace

    #Axis sizes may need changing - currently that has to be done in the script here unfortunetly - may change this later

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    ERinRiverData = pd.read_csv(Workspace + RiverPointsFileName)

    ProxCat = ERinRiverData[ERinRiverData["CatNum"] == ProxCatNo]
    DistCat = ERinRiverData[ERinRiverData["CatNum"] == DistCatNo]

    filename = 'Prox' + str(ProxCatNo) + 'Dist' + str(DistCatNo)

    # For Figure 1 - Flow distance by Elevation with color representing catchment-averaged erosion rate
    xprox = (ProxCat['FL'] / 1000)
    yprox = ProxCat['Z']
    cprox = ((ProxCat['ERvol']) / (ProxCat['FAm2'])) * 1000

    xdist = (DistCat['FL'] / 1000)
    ydist = DistCat['Z']
    cdist = ((DistCat['ERvol']) / (DistCat['FAm2'])) * 1000

    fig1, axs = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=True, figsize=[12, 6])
    cmap = cm.get_cmap('Reds')
    norm = cm.colors.Normalize()

    prox = axs[0].scatter(xprox, yprox, 25, cprox, cmap=cmap, norm=norm, edgecolors='k', linewidth=0.3)
    dist = axs[1].scatter(xdist, ydist, 25, cdist, cmap=cmap, norm=norm, edgecolors='k', linewidth=0.3)
    # ax2.invert_xaxis()

    cbar = fig1.colorbar(dist, ax=axs.ravel().tolist())
    cbar.set_label('Catchment-Averaged Erosion Rate ($mm yr^{-1}$)', fontsize=15)
    cbar.ax.tick_params(labelsize=15)
    axs[0].set_xlabel('Flow Distance to Outlet (km)', fontsize=15, weight='bold')
    axs[1].set_xlabel('Flow Distance to Outlet (km)', fontsize=15, weight='bold')
    axs[0].set_ylabel('Elevation (m)', fontsize=15, weight='bold')
    axs[0].set_xlim(0, 26)
    axs[1].set_xlim(26, 0)
    axs[0].set_ylim(0, 500)
    axs[1].set_ylim(0, 500)
    axs[0].tick_params(axis='both', labelsize=15)
    axs[1].tick_params(axis='both', labelsize=15)
    axs[0].set_xticks(range(0, 27, 2))
    axs[1].set_xticks(range(0, 27, 2))
    axs[0].text(0.5, 490, 'Proximal Flank', fontsize=15, va='top', ha='left')
    axs[1].text(0.5, 490, 'Distal Flank', fontsize=15, va='top', ha='right')

    fig1.savefig(Workspace + "fig1" + filename + ".jpg")

    # For Figure 2 - Flow distance by catchment-averaged erosion rate with color representing catchment-averaged erosion rate
    xprox = (ProxCat['FL'] / 1000)
    yprox = ((ProxCat['ERvol']) / (ProxCat['FAm2'])) * 1000
    cprox = ((ProxCat['ERvol']) / (ProxCat['FAm2'])) * 1000

    xdist = (DistCat['FL'] / 1000)
    ydist = ((DistCat['ERvol']) / (DistCat['FAm2'])) * 1000
    cdist = ((DistCat['ERvol']) / (DistCat['FAm2'])) * 1000

    fig2, axs = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=True, figsize=[12, 6])
    cmap = cm.get_cmap('Reds')
    norm = cm.colors.Normalize()

    prox = axs[0].scatter(xprox, yprox, 25, cprox, cmap=cmap, norm=norm, edgecolors='k', linewidth=0.3)
    dist = axs[1].scatter(xdist, ydist, 25, cdist, cmap=cmap, norm=norm, edgecolors='k', linewidth=0.3)
    # ax4.invert_xaxis()

    cbar = fig2.colorbar(dist, ax=axs.ravel().tolist())
    cbar.set_label('Catchment-Averaged Erosion Rate ($mm yr^{-1}$)', fontsize=15)
    cbar.ax.tick_params(labelsize=15)
    axs[0].set_xlabel('Flow Distance to Outlet (km)', fontsize=15, weight='bold')
    axs[1].set_xlabel('Flow Distance to Outlet (km)', fontsize=15, weight='bold')
    axs[0].set_ylabel('Catchment-Averaged Erosion Rate ($mm yr^{-1}$)', fontsize=15, weight='bold')
    axs[0].set_xlim(0, 26)
    axs[1].set_xlim(26, 0)
    axs[0].set_ylim(0, 0.25)
    axs[1].set_ylim(0, 0.25)
    axs[0].tick_params(axis='both', labelsize=15)
    axs[1].tick_params(axis='both', labelsize=15)
    axs[0].set_xticks(range(0, 27, 2))
    axs[1].set_xticks(range(0, 27, 2))
    axs[0].text(0.5, 0.24, 'Proximal Flank', fontsize=15, va='top', ha='left')
    axs[1].text(0.5, 0.24, 'Distal Flank', fontsize=15, va='top', ha='right')

    fig2.savefig(Workspace + "fig2" + filename + ".jpg")

    # For Figure 3 Distance to MDD by catchment-averaged erosion rate with color represenitng catchment-averaged erosion rate

    xprox = (ProxCat['NEAR_DIST'] / 1000)
    yprox = ((ProxCat['ERvol']) / (ProxCat['FAm2'])) * 1000
    cprox = ((ProxCat['ERvol']) / (ProxCat['FAm2'])) * 1000

    xdist = (DistCat['NEAR_DIST'] / 1000)
    ydist = ((DistCat['ERvol']) / (DistCat['FAm2'])) * 1000
    cdist = ((DistCat['ERvol']) / (DistCat['FAm2'])) * 1000

    fig3, axs = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=True, figsize=[12, 6])
    cmap = cm.get_cmap('Reds')
    norm = cm.colors.Normalize()

    prox = axs[0].scatter(xprox, yprox, 25, cprox, cmap=cmap, norm=norm, edgecolors='k', linewidth=0.3)
    dist = axs[1].scatter(xdist, ydist, 25, cdist, cmap=cmap, norm=norm, edgecolors='k', linewidth=0.3)
    # ax5.invert_xaxis()

    cbar = fig3.colorbar(dist, ax=axs.ravel().tolist())
    cbar.set_label('Catchment-Averaged Erosion Rate ($mm yr^{-1}$)', fontsize=15)
    cbar.ax.tick_params(labelsize=15)
    axs[0].set_xlabel('Distance from the MDD (km)', fontsize=15, weight='bold')
    axs[1].set_xlabel('Distance from the MDD (km)', fontsize=15, weight='bold')
    axs[0].set_ylabel('Catchment-Averaged Erosion Rate ($mm yr^{-1}$)', fontsize=15, weight='bold')
    axs[0].set_xlim(18, 0)
    axs[1].set_xlim(0, 18)
    axs[0].set_ylim(0, 0.25)
    axs[1].set_ylim(0, 0.25)
    axs[0].tick_params(axis='both', labelsize=15)
    axs[1].tick_params(axis='both', labelsize=15)
    axs[0].set_xticks(range(0, 19, 2))
    axs[1].set_xticks(range(0, 19, 2))
    axs[0].text(0.4, 0.24, 'Proximal Flank', fontsize=15, va='top', ha='right')
    axs[1].text(0.4, 0.24, 'Distal Flank', fontsize=15, va='top', ha='left')

    fig3.savefig(Workspace + "fig3" + filename + ".jpg")

    # For Figure 4 Distance to the MDD by elevation with color representing catchment-averaged erosion rate
    xprox = (ProxCat['NEAR_DIST'] / 1000)
    yprox = ProxCat['Z']
    cprox = ((ProxCat['ERvol']) / (ProxCat['FAm2'])) * 1000

    xdist = (DistCat['NEAR_DIST'] / 1000)
    ydist = DistCat['Z']
    cdist = ((DistCat['ERvol']) / (DistCat['FAm2'])) * 1000

    fig4, axs = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=True, figsize=[12, 6])
    cmap = cm.get_cmap('Reds')
    norm = cm.colors.Normalize()

    prox = axs[0].scatter(xprox, yprox, 25, cprox, cmap=cmap, norm=norm, edgecolors='k', linewidth=0.3)
    dist = axs[1].scatter(xdist, ydist, 25, cdist, cmap=cmap, norm=norm, edgecolors='k', linewidth=0.3)
    # ax5.invert_xaxis()

    cbar = fig4.colorbar(dist, ax=axs.ravel().tolist())
    cbar.set_label('Catchment-Averaged Erosion Rate ($mm yr^{-1}$)', fontsize=15)
    cbar.ax.tick_params(labelsize=15)
    axs[0].set_xlabel('Distance from the MDD (km)', fontsize=15, weight='bold')
    axs[1].set_xlabel('Distance from the MDD (km)', fontsize=15, weight='bold')
    axs[0].set_ylabel('Elevation (m)', fontsize=15, weight='bold')
    axs[0].set_xlim(18, 0)
    axs[1].set_xlim(0, 18)
    axs[0].set_ylim(0, 500)
    axs[1].set_ylim(0, 500)
    axs[0].tick_params(axis='both', labelsize=15)
    axs[1].tick_params(axis='both', labelsize=15)
    axs[0].set_xticks(range(0, 19, 2))
    axs[1].set_xticks(range(0, 19, 2))
    axs[0].text(0.4, 490, 'Proximal Flank', fontsize=15, va='top', ha='right')
    axs[1].text(0.4, 490, 'Distal Flank', fontsize=15, va='top', ha='left')

    fig4.savefig(Workspace + "fig4" + filename + ".jpg")

