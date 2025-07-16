/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var testingS_crop = ee.FeatureCollection("users/JY54/230801/YRDcrop2023"),
    testingS_grass = ee.FeatureCollection("users/JY54/230801/YRDgrass2023"),
    YRD = ee.FeatureCollection("projects/ee-jy54/assets/YRD"),
    S1IC = ee.ImageCollection("COPERNICUS/S1_GRD"),
    S2IC = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED"),
    L8IC = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
    Sample = ee.FeatureCollection("users/JY54/orinSample0428"),
    Blocks = ee.FeatureCollection("users/JY54/crop240429/NCPPlainBlocks"),
    geometry = /* color: #d63000 */ee.Geometry.MultiPoint();
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//var samples = corn.merge(cotton).merge(rice).merge(wheat_corn).merge(wheat_bean).merge(wheat_rice);
var startDate = "2023-01-01";
var endDate = "2024-01-01";
var year = 2023;
// var TargetROI = YRD;
var TargetROI = YRD//Blocks.filter(ee.Filter.eq('ID',20)).geometry();
Map.addLayer(TargetROI)
var Th_sos = 0.3;
var Th_eos = 0.3;
var Th_sdps = 0.7;
var Th_edps = 0.4; 
var Th_Tans = 20;


{//Function
    //CloudPT
  var rmCloudPT = function(collection,size){
    size = ee.Number(size);
    var First_image = ee.Image(collection.get(0));
    var Last_image = ee.Image(collection.get(size.subtract(1)));
    //求正向差  
    var diffImage1 = ee.List.sequence(1, size.subtract(1)).map(function(index) {
      index = ee.Number(index);
      var currentImage = ee.Image(collection.get(index));
      var previousImage = ee.Image(collection.get(index.subtract(1)));
      if (!previousImage) {return ee.Image.constant(0);}//不影响 
      return currentImage.subtract(previousImage);
    }).remove(ee.Image.constant(0));
    //求负向差 
    var diffImage2 =ee.List.sequence(1, size.subtract(1)).map(function(index) {
      index = ee.Number(index);
      var currentImage = ee.Image(collection.get(index.subtract(1)));
      var previousImage = ee.Image(collection.get(index));
      if (!currentImage) {return ee.Image.constant(0);}
      return currentImage.subtract(previousImage);
    }).remove(ee.Image.constant(0));
    //求和   
    var IC_add = function (ICL1,ICL2){
      var diff_add = ee.List.sequence(1, ICL1.size().subtract(1)).map(function(index){
        index = ee.Number(index);
        return ee.Image(ICL1.get(index)).add(ee.Image(ICL2.get(index.subtract(1))));
      });
      return diff_add;
    };
    var diff_add = IC_add(diffImage1,diffImage2);
    //去除异常点 
    var smoth =ee.ImageCollection(ee.List.sequence(1, size.subtract(2)).map(function(index){
       index = ee.Number(index).subtract(1);
       var cloudimage = ee.Image(collection.get(index.add(1))).add(ee.Image(diff_add.get(index)).divide(2));
       var replacedImage = ee.Image(collection.get(index.add(1)))
                                .where(ee.Image(diff_add.get(index)).gt(0.3), cloudimage);
      return replacedImage; 
     }));
    //还原全序列 
    var smoth_collection = smoth.merge(First_image).merge(Last_image).sort('doy');
    return smoth_collection;
  };
  //Image Function 
  {
    function s1forest_cla(img){
    var VHMask= img.select("VH").gt(-20);
    var forest_s1=img.select('VH').updateMask(VHMask).rename('forest_s1');
    return forest_s1;
  }
    var addTimeBand = function(img) {
      var DOY=ee.Date(img.get('system:time_start')).getRelative('day', 'year');
      img=img.set('DOY',DOY);
      img=img.addBands(img.metadata('DOY'));
      return img;
    };
    var addVIs = function(img){
      var nir = img.select('nir').divide(10000);
      var red = img.select('red').divide(10000);
      var blue = img.select('blue').divide(10000);
      var ndvi = img.normalizedDifference(['nir','red']).rename('NDVI');
      var lswi = img.normalizedDifference(['nir','swir']).rename('LSWI');
      img =img.addBands(ndvi).addBands(lswi);
      return img.select(['LSWI','NDVI']);
    };
    var maskS2cloud = function (image) { 
      var qa = image.select('QA60'); 
      var cloudBitMask = 1 << 10; 
      var cirrusBitMask = 1 << 11; 
      var mask = qa.bitwiseAnd(cloudBitMask).eq(0) 
                  .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
      return image.updateMask(mask);
    };
    var maskL8cloud = function(image) {
      var cloudShadowBitMask = 1 << 3;
      var cloudsBitMask = 1 << 5;
      var qa = image.select('QA_PIXEL');
      var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
          .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
      return image.updateMask(mask)
          .copyProperties(image, ["system:time_start"]);
    };
  }
  function filterByRange(collection, startDay, endDay, BandName) {
    var filteredCollection = collection.map(function(im) {
      var doy = im.select(BandName); // Select the Day of Year band
      var inRange = doy.gte(startDay).and(doy.lte(endDay)); // Check if DOY is within the range
      return im.updateMask(inRange).copyProperties(im, ['system:time_start']); // Apply the mask and copy properties
    });
    return filteredCollection;
  }
  
  function updateSeasonDates(oldDate, newDate) {
      var maskNew = newDate.mask();
      var maskedOld = oldDate.updateMask(maskNew.not());
      return maskedOld.unmask(newDate);
  }
  
}
}

{//ImageCollection
  var CollectionS2 = S2IC.filterDate(startDate, endDate)
                         .filterBounds(TargetROI)
                         .select(['B2','B3', 'B4', 'B8', 'B11','QA60'],['blue','green','red','nir','swir','QA60'])
                         .map(maskS2cloud)
                         .map(addVIs)
                         .map(addTimeBand);
  var CollectionL8 = L8IC.filterDate(startDate, endDate)
                 .filterBounds(TargetROI)
                 .select(['SR_B2','SR_B3', 'SR_B4', 'SR_B5', 'SR_B6','QA_PIXEL'],['blue','green','red','nir','swir','QA_PIXEL'])
                 .map(maskL8cloud)
                 .map(addVIs)
                 .map(addTimeBand);
  var HLS_collection =  ee.ImageCollection(CollectionS2.merge(CollectionL8)).sort('system:time_start'); 
  var S1 = S1IC.filterDate(startDate,endDate) 
               .filterBounds(TargetROI)
               .filterMetadata('instrumentMode', 'equals', 'IW')
               .select(['VV','VH'])
              // .map(function(image) {
              //     var edge = image.lt(-30.0);
              //     var maskedImage = image.mask().and(edge.not());
              //     return image.updateMask(maskedImage);
              //   })
                .map(addTimeBand)
                .sort('system:time_start');
  var s1_count = S1.select('VV').reduce(ee.Reducer.count());

var NDVImax = HLS_collection
                  .map(function(image){return image.updateMask(image.select('NDVI').lt(0.75));})
                  .select('NDVI').max().rename('NDVImax').clip(TargetROI);
var vegetation1 = ee.Image(1).updateMask(NDVImax.gt(0.3));
//forest detection
var S1_forest = S1.map(s1forest_cla);
var forests1_count = S1_forest.reduce(ee.Reducer.count());
var ff=ee.Image.cat([forests1_count,s1_count]);
function addforestfres1(img){
  var fre = img.expression("b1/b2",
    {'b1':img.select("forest_s1_count"),
  'b2':img.select("VV_count")}
  )
  .rename('forests1_fre');
  return fre;
}
var forests1_fre=addforestfres1(ff).clip(TargetROI);
forests1_fre = forests1_fre.unmask(ee.Image(0)).clip(TargetROI);
var vegetation2 = ee.Image(1).updateMask(forests1_fre.lt(0.65)).updateMask(forests1_fre.gt(0.1)).clip(TargetROI).rename("Forest");
var vegetation = vegetation2.updateMask(vegetation1.mask());
Map.addLayer(vegetation.clip(TargetROI),{min:1, max:1,palette:'#1cc219'}, 'vegetation',true);



var collection = HLS_collection;
var startDoy = ee.Date(startDate).getRelative('day','year');
var endDoy = (ee.Date(endDate).advance(-1,'day')).getRelative('day','year');
var starts = ee.List.sequence(startDoy, endDoy, 10);
var composites = starts.map(function(start) {
  var doy = start;
  var filtered = collection.filter(ee.Filter.dayOfYear(start, ee.Number(start).add(10))).mean();
  var bandLength = filtered.bandNames().length();
  var mask = ee.Algorithms.If({                   // mask must be done for time band
    condition : ee.Number(bandLength).gt(0),
    trueCase : filtered.select(0).mask(),
    falseCase : ee.Image(0).clip(TargetROI)    
  });
  return filtered.updateMask(mask)
                 .set('system:time_start',ee.Date.fromYMD(year,1,1).advance(doy,'day').millis())
                 .set('doy',doy);
});

var composites_cloud = rmCloudPT(composites,37);
//composites_cloud = rmCloudPT(composites_cloud,37); //List 问题 

Map.addLayer(composites_cloud.select('NDVI'),{},"NDVI",0);
Map.addLayer(composites_cloud.select('LSWI'),{},"LSWI",0);
print(composites_cloud)

}

{
var NDVI_10dCol = composites_cloud.select(['NDVI','DOY']);
//分区 
var maxND_Date = NDVI_10dCol.qualityMosaic('NDVI').select('DOY'); //maxNDVI 
var maxND = NDVI_10dCol.select('NDVI').max(); //maxNDVI
var Left_Col  = filterByRange(NDVI_10dCol,ee.Image(0),maxND_Date,'DOY');
var Right_Col = filterByRange(NDVI_10dCol,maxND_Date,ee.Image(365),'DOY');

// threshold
var minND_left = Left_Col.select('NDVI').min();
var amplitude_left = maxND.subtract(minND_left);
var thresh_SOS = amplitude_left.multiply(Th_sos).add(minND_left); 
var thresh_SDPS = amplitude_left.multiply(Th_sdps).add(minND_left); 

var minND_right = Right_Col.select('NDVI').min();
var amplitude_right = maxND.subtract(minND_right);
var thresh_EOS = amplitude_right.multiply(Th_eos).add(minND_right);
var thresh_EDPS = amplitude_right.multiply(Th_edps).add(minND_right); 
// feature [SOS EOS SDPS EDPS]
var col_aboveThresh_SOS = filterByRange(Left_Col,ee.Image(-999),thresh_SOS,'NDVI');
var SOS = col_aboveThresh_SOS.reduce(ee.Reducer.lastNonNull()).select('DOY_last').rename('SoS');
var col_aboveThresh_SDPS = filterByRange(Left_Col,ee.Image(-999),thresh_SDPS,'NDVI');
var SDPS = col_aboveThresh_SDPS.reduce(ee.Reducer.lastNonNull()).select('DOY_last').rename('SDPS');

var col_aboveThresh_EOS = filterByRange(Right_Col,ee.Image(-999),thresh_EOS,'NDVI');
var EOS = col_aboveThresh_EOS.reduce(ee.Reducer.firstNonNull()).select('DOY_first').rename('EoS');
var col_aboveThresh_EDPS = filterByRange(Right_Col,ee.Image(-999),thresh_EDPS,'NDVI');
var EDPS = col_aboveThresh_EDPS.reduce(ee.Reducer.firstNonNull()).select('DOY_first').rename('EDPS');

//双峰 

var OtherSeason_Col = NDVI_10dCol.map(function(im) {
  var condition = maxND_Date.select('DOY').gt(160);
  var maskDOYGT160 = im.select('DOY').lte(SOS); // 如果DOY大于183
  var maskDOYLTE160 = im.select('DOY').gte(EOS); // 如果DOY小于等于183
  var mask = condition.where(condition, maskDOYGT160).where(condition.not(), maskDOYLTE160);
  return im.updateMask(mask);
});
var OtherSeason_maxND = OtherSeason_Col.qualityMosaic('NDVI');
var SecPeak = OtherSeason_Col.qualityMosaic('NDVI').select('NDVI').rename('SecPeak');//.select('DOY');//Feature1

{
  var maxND_Date2 = OtherSeason_Col.qualityMosaic('NDVI').select('DOY'); //maxNDVI 
  var maxND2 = OtherSeason_Col.select('NDVI').max(); //maxNDVI
  var mask = maxND_Date.lte(160).updateMask(SecPeak.gte(0.4));
  OtherSeason_Col = OtherSeason_Col.map(function(im){ return im.updateMask(mask); });
  
  var Left_Col2 = filterByRange(OtherSeason_Col, ee.Image(0), maxND_Date2, 'DOY');
  var Right_Col2 = filterByRange(OtherSeason_Col, maxND_Date2, ee.Image(365), 'DOY');

  
  var minND_left2 = Left_Col2.select('NDVI').min();
  var amplitude_left2 = maxND2.subtract(minND_left2);
  var thresh_SOS2 = amplitude_left2.multiply(Th_sos).add(minND_left2); //th = 0.3
  var thresh_SDPS2 = amplitude_left2.multiply(Th_sdps).add(minND_left2);
  
  var minND_right2 = Right_Col2.select('NDVI').min();
  var amplitude_right2 = maxND2.subtract(minND_right2);
  var thresh_EOS2 = amplitude_right2.multiply(Th_eos).add(minND_right2);
  var thresh_EDPS2 = amplitude_right2.multiply(Th_sdps).add(minND_right2);
  

  var col_aboveThresh_SOS2 = filterByRange(Left_Col2,ee.Image(-999),thresh_SOS2,'NDVI');
  var col_aboveThresh_SDPS2 = filterByRange(Left_Col2,ee.Image(-999),thresh_SDPS2,'NDVI');
  var SOS2 = col_aboveThresh_SOS2.reduce(ee.Reducer.lastNonNull()).select('DOY_last').rename('SoS');
  var SDPS2 = col_aboveThresh_SDPS2.reduce(ee.Reducer.lastNonNull()).select('DOY_last').rename('SDPS');
 
  var col_aboveThresh_EOS2 = filterByRange(Right_Col2,ee.Image(-999),thresh_EOS2,'NDVI');
  var col_aboveThresh_EDPS2 = filterByRange(Right_Col2,ee.Image(-999),thresh_EDPS2,'NDVI'); 
  var EOS2 = col_aboveThresh_EOS2.reduce(ee.Reducer.firstNonNull()).select('DOY_first').rename('EoS');//.metadata('system:time_start','date1');
  var EDPS2 = col_aboveThresh_EDPS2.reduce(ee.Reducer.firstNonNull()).select('DOY_first').rename('EDPS');
}

var SOS = updateSeasonDates(SOS,SOS2);
var EOS = updateSeasonDates(EOS,EOS2);
var SDPS = updateSeasonDates(SDPS,SDPS2);
var EDPS = updateSeasonDates(EDPS,EDPS2);
var thresh_SOS = updateSeasonDates(thresh_SOS,thresh_SOS2);
var thresh_SDPS = updateSeasonDates(thresh_SDPS,thresh_SDPS2);

var col_MainSeason = filterByRange(composites_cloud,SOS,EOS,'DOY').select(['NDVI','LSWI']); 
print(SOS)

var VV_GSmean = filterByRange(S1,SDPS,EDPS,'DOY').select('VV').mean().unmask(-10).rename('VV_GSmean');

var LSWI_Transmean = composites_cloud.map(function(im){
  var out = im.select('DOY').gte(SOS.subtract(ee.Image(Th_Tans))).and(im.select('DOY').lte(SOS));
  var LSWI_NDVIdiff = im.select('NDVI').subtract(im.select('LSWI')).rename('TransF');
  return LSWI_NDVIdiff.updateMask(out);
}).select('TransF');


var GSL = EOS.subtract(SOS).rename('GSL');
var cropmask =vegetation.updateMask(GSL.lt(160))
Map.addLayer(cropmask.clip(TargetROI),{min:0, max:0,palette:'#1cc219'}, 'cropmask',true);
var GUS = (thresh_SDPS.subtract(thresh_SOS)).divide(SDPS.subtract(SOS)).rename('GUS');
}
// var features = VV_GSmean.addBands(LSWI_Transmean.sum()).addBands(SecPeak).addBands(GUS);
var features = VV_GSmean.addBands(LSWI_Transmean.sum()).addBands(SecPeak).addBands(GUS).addBands(SOS).addBands(EOS).addBands(GSL).addBands(SDPS);
Map.addLayer(features,{},'Features')
print(Sample.limit(10))

// var randomPoints = ee.FeatureCollection.randomPoints({
//   region:YRD,//vegetation.clip(TargetROI).geometry(),
//   points:1000,
//   seed:20
// });

// Map.addLayer(randomPoints, {color: 'red'}, 'Random Points');

var CropProperties = ["GUS","SecPeak","TransF","VV_GSmean"]
// var RPoints = features.sampleRegions({
//     collection: randomPoints,
//     // properties: inputProperties,
//     scale: 10,
//     tileScale:4,
//     geometries:true
// })
// print(RPoints.limit(100))
// var classified0 = RPoints.classify(rf_classifier);
// print(classified0)
// var classProbabilities = RPoints.classify(rf_classifier.setOutputMode('MULTIPROBABILITY'));
// print(classProbabilities)

var rf_classifier= ee.Classifier.smileRandomForest(150).train(Sample, 'class',CropProperties,0.5,20);
print(rf_classifier.explain())


var classified1 = features.classify(rf_classifier);
var PROBABILITY = features.classify(rf_classifier.setOutputMode('MULTIPROBABILITY'));
var bands = ['corn','cotton','rice','wheat_corn','wheat_bean','wheat_rice']
var cla_prob = PROBABILITY.arrayFlatten([bands])
var cla_probIC = ee.ImageCollection.fromImages(bands.map(function(n){return cla_prob.select([n]).rename('prob')}));
print(cla_probIC)
// var cla_confidentIC = cla_probIC.map(function(image){
//   return image.updateMask(image.gt(0.7))
// }).mean().clip(TargetROI).updateMask(vegetation)
// print(cla_confidentIC)
// Map.addLayer(cla_confidentIC)
var cla_confidentIC = cla_probIC.max().updateMask(cla_probIC.max().gt(0.7)).clip(TargetROI).updateMask(cropmask)
var cla_sampleimg = classified1.updateMask(cla_confidentIC)
print(cla_sampleimg)
var samples =cla_sampleimg.stratifiedSample({
  numPoints:100,
  classBand:'classification',
  region:TargetROI,
  scale:100,
  geometries:true,
  tileScale:16
  
})
print(samples)

var newSample = features.sampleRegions({collection: samples, properties: ['classification'], scale: 30,tileScale:16,geometries:true});
print(newSample,'newSample')
Export.table.toDrive(newSample)
print(cla_confidentIC)
Map.addLayer(cla_probIC.max(),{},'cla_probIC')
// Map.addLayer(cla_confidentIC.mean().clip(TargetROI).updateMask(vegetation))
var palette = ['1f77b4', 'ff7f0e','2ca02c',   'd62728',   '9467bd',   '8c564b'];
Map.addLayer(classified1.updateMask(cropmask).clip(TargetROI),{max:6,min:1,palette:palette},'classified1')
Map.addLayer(cla_prob.updateMask(cropmask).clip(TargetROI))

var RegionProperties = ["VV_GSmean","TransF","SecPeak","GUS","SoS","EoS","GSL","SDPS"]

