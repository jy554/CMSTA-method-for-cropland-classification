/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var features = ee.Image("users/JY54/crop240429/features"),
    newSample = ee.FeatureCollection("users/JY54/crop240429/targetSample"),
    Sample = ee.FeatureCollection("users/JY54/crop240429/SourceSample"),
    NCPBlocks = ee.FeatureCollection("users/JY54/crop240429/NCPPlainBlocks"),
    YRD = ee.FeatureCollection("projects/ee-jy54/assets/YRD");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// Export.table.toDrive(newSample)
// // Export.table.toDrive(Sample)
// Export.image.toDrive({image:features.float(),maxPixels:10000000000000})

var id = "YRD0805";
var fpath = "users/JY54/NCPFeatures/NCPFeature"+"YRD";
var Ppath = "users/JY54/NCPSamples/NCPSample"+id;
features = ee.Image(fpath);
newSample = ee.FeatureCollection(Ppath);
var ROI = YRD;//NCPBlocks.filter(ee.Filter.eq('ID',id)).geometry();

//TrAdaboost function
{
  function TrAdaboost(FeatureImage,SourceSample,TargetSample,numIterations){
    var samplenum = ee.Number(1e5);
    //初始化权重,beta 
    var sourceWeights = ee.List.repeat(ee.Number(0.5).divide(SourceSample.size()), SourceSample.size());
    var targetWeights = ee.List.repeat(ee.Number(0.5).divide(TargetSample.size()), TargetSample.size());
    SourceSample = CalnumByweights(addWeight(addID(SourceSample),sourceWeights));
    TargetSample = CalnumByweights(addWeight(addID(TargetSample),targetWeights));
    var beta = ee.Number.expression('1/(1+(2*log(n)/N)**0.5)',{n:SourceSample.size(),N:numIterations});
    //迭代
    for(var i = 0;i < numIterations;i++){
      //利用调整样本数量的方式调整样本权值
      TargetSample = CopyFeatureByweight(TargetSample);
      SourceSample = CopyFeatureByweight(SourceSample);
      var Sample = ee.FeatureCollection([TargetSample,SourceSample]).flatten();
      var trainedClassifier = ee.Classifier.minimumDistance().train(Sample,'class',FeatureImage.bandNames(),1,0); 
      Sample = Sample.distinct('ID');
      TargetSample = TargetSample.distinct('ID');
      var trainedclassified_error = TargetSample.classify(trainedClassifier);
      var trainedclassified = Sample.classify(trainedClassifier);
      var error = calculate_error_rate(trainedclassified_error,TargetSample);
      //计算适应性权重                     
      var beta_T = ee.Number(error).divide(ee.Number(1).subtract(error)); 
      //更新权重 
      SourceSample = trainedclassified.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);
      TargetSample = trainedclassified.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
      SourceSample = CalnumByweights(UpdateWeights(SourceSample,beta,beta_T),samplenum);
      TargetSample = CalnumByweights(UpdateWeights(TargetSample,beta,beta_T),samplenum);//.distinct('ID')
    }  
    TargetSample = CopyFeatureByweight(TargetSample);
    SourceSample = CopyFeatureByweight(SourceSample);
    var SampleSet = ee.FeatureCollection([SourceSample,TargetSample]);
    var classifier = ee.Classifier.smileRandomForest(30).train(SampleSet,'class',FeatureImage.bandNames(),1,0); 
    var classified = FeatureImage.classify(classifier);
    return classified;
  }
  //为样本添加ID属性 
  function addID(featurecollction){  
    var fc = featurecollction.map(function(feature) {
      var addSeq = function(feature, list) {
        var index = ee.List(list).indexOf(feature);
        return feature.set('ID', index);
      };
      return addSeq(feature, featurecollction.toList(featurecollction.size()));
    });
    return fc;
  }
  //为样本添加初始权重属性
  function addWeight(featurecollection,weightProperty) {
    var fc = featurecollection.map(function(feature){
       var index = ee.Feature(feature).getNumber('ID');
       var newValue = weightProperty.get(index);
       return feature.set('weight', newValue);
    });
    return fc;
  }
  //计算错误率 
  function calculate_error_rate(featurecollection,TargetSample){
    var fc = featurecollection.map(function(feature){
        // var TargetSample_cal = TargetSample.map(function(feature){
        //   var weight_cal = feature.getNumber('num').multiply(feature.getNumber('weight'));
        //   return feature.set('weight_cal',weight_cal);
        //   });
        var weight_sum = TargetSample.aggregate_sum('weight');
      var error = ee.Number.expression("weight_i/weight_sum*(y_T-y_P)",{
        weight_i:feature.getNumber('weight'),
        weight_sum:weight_sum.multiply(0.5),
        y_T:ee.Number(0.5),//feature.getNumber('class'),
        y_P:feature.getArray('classification').get([ee.Number.parse(feature.get('class')).subtract(1)]).min(0.5)//feature.getNumber('classification') 
       });
      return feature.set('error',error);
    });
    // var fc_cal = fc.map(function(feature){
    //   var error_cal = feature.getNumber('num').multiply(feature.getNumber('error'));
    //   return feature.set('error_cal',error_cal);
    // });
    return fc.aggregate_sum('error');
  }
  //调整样本权重 （ source 样本权重 正确样本权重变大;target 样本权重 错误样本权重变大 ）
  function UpdateWeights(featurecollection,beta,bata_T){
    var fc =featurecollection.map(function(feature){
      var expression = ee.Algorithms.If({
        condition:ee.Algorithms.IsEqual(feature.get('type'),'1'),//1:source
        trueCase:"(beta**(abs(y_T-y_P)*2))*weight",//Source
        falseCase:"(beta**(-abs(y_T-y_P)*2))*weight"
      });
      var beta_UW = ee.Algorithms.If({
        condition: ee.Algorithms.IsEqual(feature.get('type'),'1'),
        trueCase: beta,
        falseCase: bata_T
      });
     var Newweight=ee.Number.expression(expression,{
         beta: beta_UW,
         weight:feature.get('weight'),
         y_T:ee.Number(0.5),
         y_P:feature.getArray('classification').get([ee.Number.parse(feature.get('class')).subtract(1)]).min(0.5)
        });
      return feature.set('weight',Newweight); 
     });
     return fc;
  }
  //根据权重计算样本数量  
  function CalnumByweights(featurecollection,samplenum,weight_sum){
     var sum = ee.Number.parse(weight_sum);
     var fc = featurecollection.map(function(feature){
       var weight = feature.getNumber('weight');
       var num = ee.Algorithms.If({
        condition: ee.Algorithms.IsEqual(feature.get('type'),'1'),
        trueCase: weight.divide(sum).multiply(samplenum).round(),//.max(ee.Number(1)),//round(),ceil(),floor()
        falseCase: weight.divide(sum).multiply(samplenum).ceil()//round(),ceil()
        });
       return feature.set('num',num);
     });
     return fc; 
   }
  //根据样本数量制作数据集 
  function CopyFeatureByweight(featurecollection){
    var flattened = featurecollection.toList(featurecollection.size()).map(function(feature){
     var fe=ee.Feature(feature);
     return ee.List.repeat(fe, fe.get('num'));
    }).flatten();
    var fc = ee.FeatureCollection(flattened);
    return fc;
  }
}
var transformFeatures = function(featureCollection,Vector,vectorID) {
  // 获取原始特征的属性值
  var newfc = featureCollection.map(function(feature){
    var originalProperties = feature.toDictionary();
    var newProperties = originalProperties.map(function(key, value) {
      var newValue = ee.Algorithms.If(
        ee.String(key).match('(SoS|SDPS)'), 
        ee.Number(value).multiply(Vector[0]),
        ee.Algorithms.If(
          ee.String(key).match('GSL|EoS'), 
          ee.Number(value).multiply(Vector[1]),
          ee.Algorithms.If(
            ee.String(key).match('GUS'),
            ee.Number(value).multiply(Vector[2]),
            value)));
      return newValue;
    });
    return ee.Feature(feature.geometry(), newProperties).set('vectorType',vectorID+1);
  });
  return newfc;
};
function renameProperty(feature) {
  var oldValue = feature.get('classification');
  return feature.set({'class': oldValue}).set({'type':'2'})
}


var changeVectors = [
  [0.9, 0.9, 1.1], // 特征向量 1
  [0.9, 1.1, 0.9], // 特征向量 2
  [1.1, 0.9, 1.1], // 特征向量 3
  [1.1, 1.1, 0.9]  // 特征向量 4
  // 在这里添加更多特征向量...
  ];


var samples =  Sample.filter(ee.Filter.eq('class',1)).limit(1000)
              .merge( Sample.filter(ee.Filter.eq('class',2)).limit(1000))
              .merge( Sample.filter(ee.Filter.eq('class',3)).limit(1000))
              .merge( Sample.filter(ee.Filter.eq('class',4)).limit(1000))
              .merge( Sample.filter(ee.Filter.eq('class',5)).limit(1000))
              .merge( Sample.filter(ee.Filter.eq('class',6)).limit(1000))
var change1sample = transformFeatures(samples,changeVectors[0],0)
var change2sample = transformFeatures(samples,changeVectors[1],1)
var change3sample = transformFeatures(samples,changeVectors[2],2)
var change4sample = transformFeatures(samples,changeVectors[3],3)
print(change4sample.limit(100))
samples = samples.map(function(feature){return feature.set('vectorType',0)})
var SourceSamples = ee.FeatureCollection([samples,change1sample,change2sample,change3sample,change4sample]).flatten()//,change1sample,change2sample,change3sample,change4sample
  .map(function(feature){return feature.set('type','1')});
print(SourceSamples.size())
newSample = newSample.map(renameProperty)


{
  var bandnames = ['GUS','SoS','EoS','GSL','SDPS']//features.bandNames();//

  var numIterations = 6;
  var samplenum = ee.Number(5e4);
  var SourceSample = SourceSamples
  var TargetSample = newSample
  var FeatureImage = features;
  
  var sourceWeights = ee.List.repeat(ee.Number(0.5).divide(SourceSample.size()), SourceSample.size());
  var targetWeights = ee.List.repeat(ee.Number(0.5).divide(TargetSample.size()), TargetSample.size());
  
  SourceSample = addWeight(addID(SourceSample),sourceWeights);
  TargetSample = addWeight(addID(TargetSample),targetWeights);
  
  var weight_sum = ee.FeatureCollection([TargetSample,SourceSample]).flatten().aggregate_sum('weight');
 
  SourceSample = CalnumByweights(SourceSample,samplenum,weight_sum);
  TargetSample = CalnumByweights(TargetSample,samplenum,weight_sum);
  var beta = ee.Number.expression('1/(1+(2*log(n/N))**0.5)',{n:SourceSample.size(),N:numIterations});
  

  //-------------------------------------迭代-----------------------------------------------------
  //for(var i = 0;i < numIterations;i++){
  {  
    //111111111111111111111
    {//利用调整样本数量的方式调整样本权值
    var TargetSample1 = CopyFeatureByweight(TargetSample);
    var SourceSample1 = CopyFeatureByweight(SourceSample);
    
    var Sample1 = ee.FeatureCollection([TargetSample1,SourceSample1]).flatten();
    var trainedClassifier1 = ee.Classifier.smileCart(10).train(Sample1,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    var trainedclassified_error1 = TargetSample1.distinct('ID').classify(trainedClassifier1);
    var trainedclassified1 = Sample1.distinct('ID').classify(trainedClassifier1);
    var error1 = calculate_error_rate(trainedclassified_error1,TargetSample1.distinct('ID'));
    print(error1,'error1')
    //计算适应性权重
    error1 = error1.min(ee.Number(0.499));
    var beta_T1 = ee.Number(error1).divide(ee.Number(1).subtract(error1));
    //更新权重 
    SourceSample1 = trainedclassified1.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);
    TargetSample1 = trainedclassified1.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample1 = UpdateWeights(SourceSample1,beta,beta_T1);
    TargetSample1 = UpdateWeights(TargetSample1,beta,beta_T1);
    var weight_sum1 = ee.FeatureCollection([TargetSample1,SourceSample1]).flatten().aggregate_sum('weight');
    SourceSample1 = CalnumByweights(SourceSample1,samplenum,weight_sum1);
    TargetSample1 = CalnumByweights(TargetSample1,samplenum,weight_sum1);//.distinct('ID')
    print(weight_sum1)
    }
    //22222222222222222222
    {var TargetSample2 = CopyFeatureByweight(TargetSample1);
    var SourceSample2 = CopyFeatureByweight(SourceSample1);
    SourceSample1 = null; TargetSample1 = null;Sample1 = null;
    print(SourceSample2.size(),'SourceSample2')
    var Sample2 = ee.FeatureCollection([TargetSample2,SourceSample2]).flatten();
    var trainedClassifier2 = ee.Classifier.smileCart(20).train(Sample2,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    var trainedclassified_error2 = TargetSample2.distinct('ID').classify(trainedClassifier2);
    var trainedclassified2 = Sample2.distinct('ID').classify(trainedClassifier2);
    var error2 = calculate_error_rate(trainedclassified_error2,TargetSample2.distinct('ID'));
    print(error2,"error2");
    //计算适应性权重    
    error2 = error2.min(ee.Number(0.499));
    var beta_T2 = ee.Number(error2).divide(ee.Number(1).subtract(error2)); 
    //更新权重 
    SourceSample2 = trainedclassified2.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);
    TargetSample2 = trainedclassified2.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample2 = UpdateWeights(SourceSample2,beta,beta_T2);
    TargetSample2 = UpdateWeights(TargetSample2,beta,beta_T2);
    var weight_sum2 = ee.FeatureCollection([TargetSample2,SourceSample2]).flatten().aggregate_sum('weight');
    print(weight_sum2)
    SourceSample2 = CalnumByweights(SourceSample2,samplenum,weight_sum2);
    TargetSample2 = CalnumByweights(TargetSample2,samplenum,weight_sum2);//.distinct('ID')
    }
    //333333333333333333333
    {var TargetSample3 = CopyFeatureByweight(TargetSample2);
    var SourceSample3 = CopyFeatureByweight(SourceSample2);
    TargetSample2 = null; TargetSample2 = null;Sample2 = null;
    print(SourceSample3.size(),'SourceSample3')
    print(TargetSample3.size(),'TargetSample3')
    var Sample3 = ee.FeatureCollection([TargetSample3,SourceSample3]).flatten();
    var trainedClassifier3 = ee.Classifier.smileCart(30).train(Sample3,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    var trainedclassified_error3 = TargetSample3.distinct('ID').classify(trainedClassifier3);
    var trainedclassified3 = Sample3.distinct('ID').classify(trainedClassifier3);///////////////limit
    var error3 = calculate_error_rate(trainedclassified_error3,TargetSample3.distinct('ID'));
    print(error3,"error3");
    //计算适应性权重                     
    error3 = error3.min(ee.Number(0.499));
    var beta_T3 = ee.Number(error3).divide(ee.Number(1).subtract(error3)); 
    //更新权重 
    SourceSample3 = trainedclassified3.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample3 = trainedclassified3.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample3 = UpdateWeights(SourceSample3,beta,beta_T3);
    TargetSample3 = UpdateWeights(TargetSample3,beta,beta_T3);
    var weight_sum3 = ee.FeatureCollection([TargetSample3,SourceSample3]).flatten().aggregate_sum('weight');
    SourceSample3 = CalnumByweights(SourceSample3,samplenum,weight_sum3);
    TargetSample3 = CalnumByweights(TargetSample3,samplenum,weight_sum3);//.distinct('ID')
    }
    //444444444444444444444
    {var TargetSample4 = CopyFeatureByweight(TargetSample3);
    var SourceSample4 = CopyFeatureByweight(SourceSample3);
    SourceSample3 = null; TargetSample3 = null;Sample3 = null;
    print(SourceSample4.size(),'SourceSample4');
    var Sample4 = ee.FeatureCollection([TargetSample4,SourceSample4]).flatten();
    var trainedClassifier4 = ee.Classifier.smileCart(40).train(Sample4,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    Sample4 = Sample4.distinct('ID');
    TargetSample4 = TargetSample4.distinct('ID');
    var trainedclassified_error4 = TargetSample4.classify(trainedClassifier4);
    var trainedclassified4 = Sample4.classify(trainedClassifier4);///////////////limit
    var error4 = calculate_error_rate(trainedclassified_error4,TargetSample4.distinct('ID'));
    print(error4,"error4");
    //计算适应性权重   
    error4 = error4.min(ee.Number(0.499));
    var beta_T4 = ee.Number(error4).divide(ee.Number(1).subtract(error4)); 
    //更新权重 
    SourceSample4 = trainedclassified4.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample4 = trainedclassified4.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample4 = UpdateWeights(SourceSample4,beta,beta_T4);
    TargetSample4 = UpdateWeights(TargetSample4,beta,beta_T4);
    var weight_sum4 = ee.FeatureCollection([TargetSample4,SourceSample4]).flatten().aggregate_sum('weight');
    SourceSample4 = CalnumByweights(SourceSample4,samplenum,weight_sum4);
    TargetSample4 = CalnumByweights(TargetSample4,samplenum,weight_sum4);//.distinct('ID')
    }
    //5555555555555555555555
    {var TargetSample5 = CopyFeatureByweight(TargetSample4);
    var SourceSample5 = CopyFeatureByweight(SourceSample4);
    SourceSample4 = null; TargetSample4 = null;Sample4 = null;
    print(SourceSample5.size(),'SourceSample5')
    var Sample5 = ee.FeatureCollection([TargetSample5,SourceSample5]).flatten();
    var trainedClassifier5 = ee.Classifier.smileCart(50).train(Sample5,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    Sample5 = Sample5.distinct('ID');
    TargetSample5 = TargetSample5.distinct('ID');
    var trainedclassified_error5 = TargetSample5.classify(trainedClassifier5);
    var trainedclassified5 = Sample5.classify(trainedClassifier5);///////////////limit
    var error5 = calculate_error_rate(trainedclassified_error5,TargetSample5.distinct('ID'));
    print(error5,"error5");
    //计算适应性权重   
    error5 = error5.min(ee.Number(0.499));
    var beta_T5 = ee.Number(error5).divide(ee.Number(1).subtract(error5)); 
    //更新权重 
    SourceSample5 = trainedclassified5.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample5 = trainedclassified5.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample5 = UpdateWeights(SourceSample5,beta,beta_T5);
    TargetSample5 = UpdateWeights(TargetSample5,beta,beta_T5);
    var weight_sum5 = ee.FeatureCollection([TargetSample5,SourceSample5]).flatten().aggregate_sum('weight');
    SourceSample5 = CalnumByweights(SourceSample5,samplenum,weight_sum5);
    TargetSample5 = CalnumByweights(TargetSample5,samplenum,weight_sum5);//.distinct('ID')
    }
    //666666666666666666666666
    {var TargetSample6 = CopyFeatureByweight(TargetSample5);
    var SourceSample6 = CopyFeatureByweight(SourceSample5);
    SourceSample5 = null; TargetSample5 = null;Sample5 = null;
    var Sample6 = ee.FeatureCollection([TargetSample6,SourceSample6]).flatten();
    var trainedClassifier6 = ee.Classifier.smileCart(60).train(Sample6,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    Sample6 = Sample6.distinct('ID');
    TargetSample6 = TargetSample6.distinct('ID');
    SourceSample6 = SourceSample6.distinct('ID');
    var trainedclassified_error6 = TargetSample6.classify(trainedClassifier6);
    var trainedclassified6 = Sample6.classify(trainedClassifier6);///////////////limit
    var error6 = calculate_error_rate(trainedclassified_error6,TargetSample6.distinct('ID'));
    print(error6,"error6");
    //计算适应性权重   
    error6 = error6.min(ee.Number(0.499));
    var beta_T6 = ee.Number(error6).divide(ee.Number(1).subtract(error6)); 
    //更新权重 
    SourceSample6 = trainedclassified6.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample6 = trainedclassified6.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample6 = UpdateWeights(SourceSample6,beta,beta_T6);
    TargetSample6 = UpdateWeights(TargetSample6,beta,beta_T6);
    var weight_sum6 = ee.FeatureCollection([TargetSample6,SourceSample6]).flatten().aggregate_sum('weight');

    SourceSample6 = CalnumByweights(SourceSample6,samplenum,weight_sum6);
    TargetSample6 = CalnumByweights(TargetSample6,samplenum,weight_sum6);//.distinct('ID')
    }
    //77777777777777777777777
    {var TargetSample7 = CopyFeatureByweight(TargetSample6);
    var SourceSample7 = CopyFeatureByweight(SourceSample6);
    SourceSample6 = null; TargetSample6 = null;Sample6 = null;
    var Sample7 = ee.FeatureCollection([TargetSample7,SourceSample7]).flatten();
    print(SourceSample7.size())
    print(SourceSample7)
    var trainedClassifier7 = ee.Classifier.smileCart(70).train(Sample7,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    Sample7 = Sample7.distinct('ID');
    TargetSample7 = TargetSample7.distinct('ID');
    var trainedclassified_error7 = TargetSample7.classify(trainedClassifier7);
    var trainedclassified7 = Sample7.classify(trainedClassifier7);///////////////limit
    var error7 = calculate_error_rate(trainedclassified_error7,TargetSample7.distinct('ID'));
    print(error7,"error7");
    //计算适应性权重   
    error7 = error7.min(ee.Number(0.499));
    var beta_T7 = ee.Number(error7).divide(ee.Number(1).subtract(error7)); 
    //更新权重 
    SourceSample7 = trainedclassified7.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample7 = trainedclassified7.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample7 = UpdateWeights(SourceSample7,beta,beta_T7);
    TargetSample7 = UpdateWeights(TargetSample7,beta,beta_T7);
    var weight_sum7 = ee.FeatureCollection([TargetSample7,SourceSample7]).flatten().aggregate_sum('weight');

    SourceSample7 = CalnumByweights(SourceSample7,samplenum,weight_sum7);
    TargetSample7 = CalnumByweights(TargetSample7,samplenum,weight_sum7);//.distinct('ID')
    }
    //888888888888888888888888
    {var TargetSample8 = CopyFeatureByweight(TargetSample7);
    var SourceSample8 = CopyFeatureByweight(SourceSample7);
    SourceSample7 = null; TargetSample7 = null;
    var Sample8 = ee.FeatureCollection([TargetSample8,SourceSample8]).flatten();
    var trainedClassifier8 = ee.Classifier.smileCart(80).train(Sample8,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    Sample8 = Sample8.distinct('ID');
    TargetSample8 = TargetSample8.distinct('ID');
    var trainedclassified_error8 = TargetSample8.classify(trainedClassifier8);
    var trainedclassified8 = Sample8.classify(trainedClassifier8);///////////////limit
    var error8 = calculate_error_rate(trainedclassified_error8,TargetSample8.distinct('ID'));
    print(error8,"error8");
    //计算适应性权重   
    error8 = error8.min(ee.Number(0.499));
    var beta_T8 = ee.Number(error8).divide(ee.Number(1).subtract(error8)); 
    //更新权重 
    SourceSample8 = trainedclassified8.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample8 = trainedclassified8.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample8 = UpdateWeights(SourceSample8,beta,beta_T8);
    TargetSample8 = UpdateWeights(TargetSample8,beta,beta_T8);
    var weight_sum8 = ee.FeatureCollection([TargetSample8,SourceSample8]).flatten().aggregate_sum('weight');

    SourceSample8 = CalnumByweights(SourceSample8,samplenum,weight_sum8);
    TargetSample8 = CalnumByweights(TargetSample8,samplenum,weight_sum8);//.distinct('ID')
    }
    //9999999999999999999999999
    {var TargetSample9 = CopyFeatureByweight(TargetSample8);
    var SourceSample9 = CopyFeatureByweight(SourceSample8);
    SourceSample8 = null; TargetSample8 = null;Sample8 = null;
    var Sample9 = ee.FeatureCollection([TargetSample9,SourceSample9]).flatten();
    var trainedClassifier9 = ee.Classifier.smileCart(90).train(Sample9,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    Sample9 = Sample9.distinct('ID');
    TargetSample9 = TargetSample9.distinct('ID');
    var trainedclassified_error9 = TargetSample9.classify(trainedClassifier9);
    var trainedclassified9 = Sample9.classify(trainedClassifier9);
    var error9 = calculate_error_rate(trainedclassified_error9,TargetSample9.distinct('ID'));
    print(error9,"error9");
    //计算适应性权重   
    error9 = error9.min(ee.Number(0.499));
    var beta_T9 = ee.Number(error9).divide(ee.Number(1).subtract(error9)); 
    //更新权重 
    SourceSample9 = trainedclassified9.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample9 = trainedclassified9.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample9 = UpdateWeights(SourceSample9,beta,beta_T9);
    TargetSample9 = UpdateWeights(TargetSample9,beta,beta_T9);
    var weight_sum9 = ee.FeatureCollection([TargetSample9,SourceSample9]).flatten().aggregate_sum('weight');

    SourceSample9 = CalnumByweights(SourceSample9,samplenum,weight_sum9);
    TargetSample9 = CalnumByweights(TargetSample9,samplenum,weight_sum9);//.distinct('ID')
    }
    //10
    {var TargetSample10 = CopyFeatureByweight(TargetSample9);
    var SourceSample10 = CopyFeatureByweight(SourceSample9);
    SourceSample9 = null; TargetSample9 = null;Sample9 = null;
    var Sample10 = ee.FeatureCollection([TargetSample10,SourceSample10]).flatten();
    var trainedClassifier10 = ee.Classifier.smileCart(100).train(Sample10,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    Sample10 = Sample10.distinct('ID');
    TargetSample10 = TargetSample10.distinct('ID');
    var trainedclassified_error10 = TargetSample10.classify(trainedClassifier10);
    var trainedclassified10 = Sample10.classify(trainedClassifier10);
    var error10 = calculate_error_rate(trainedclassified_error10,TargetSample10.distinct('ID'));
    print(error10,"error10");
    //计算适应性权重   
    error10 = error10.min(ee.Number(0.499));
    var beta_T10 = ee.Number(error10).divide(ee.Number(1).subtract(error10)); 
    //更新权重 
    SourceSample10 = trainedclassified10.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample10 = trainedclassified10.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample10 = UpdateWeights(SourceSample10,beta,beta_T10);
    TargetSample10 = UpdateWeights(TargetSample10,beta,beta_T10);
    var weight_sum10 = ee.FeatureCollection([TargetSample10,SourceSample10]).flatten().aggregate_sum('weight');

    SourceSample10 = CalnumByweights(SourceSample10,samplenum,weight_sum10);
    TargetSample10 = CalnumByweights(TargetSample10,samplenum,weight_sum10);//.distinct('ID')
    }
    //11
    {var TargetSample11 = CopyFeatureByweight(TargetSample10);
    var SourceSample11 = CopyFeatureByweight(SourceSample10);
    SourceSample10 = null; TargetSample10 = null;Sample10 = null;
    var Sample11 = ee.FeatureCollection([TargetSample11,SourceSample11]).flatten();
    var trainedClassifier11 = ee.Classifier.smileCart(100).train(Sample11,'class',bandnames,1,0).setOutputMode('MULTIPROBABILITY'); 
    Sample11 = Sample11.distinct('ID');
    TargetSample11 = TargetSample11.distinct('ID');
    var trainedclassified_error11 = TargetSample11.classify(trainedClassifier11);
    var trainedclassified11 = Sample11.classify(trainedClassifier11);
    var error11 = calculate_error_rate(trainedclassified_error11,TargetSample11.distinct('ID'));
    print(error11,"error11");
    //计算适应性权重   
    error11 = error11.min(ee.Number(0.499));
    var beta_T11 = ee.Number(error11).divide(ee.Number(1).subtract(error11)); 
    //更新权重 
    SourceSample11 = trainedclassified11.filter(ee.Filter.eq('type','1'));//filterMetadata('type',"equals",2);/////////limit
    TargetSample11 = trainedclassified11.filter(ee.Filter.eq('type','2'));//.filterMetadata('type',"equals",1);
    SourceSample11 = UpdateWeights(SourceSample11,beta,beta_T11);
    TargetSample11 = UpdateWeights(TargetSample11,beta,beta_T11);
    var weight_sum11 = ee.FeatureCollection([TargetSample11,SourceSample11]).flatten().aggregate_sum('weight');

    SourceSample11 = CalnumByweights(SourceSample11,samplenum,weight_sum11);
    TargetSample11 = CalnumByweights(TargetSample11,samplenum,weight_sum11);//.distinct('ID')
    }
  }  
  //}  
  print(SourceSample11.distinct('ID').size())

  
  var TargetSample_cls = CopyFeatureByweight(SourceSample11);
  var SourceSample_cls = CopyFeatureByweight(SourceSample11);//SourceSample6 //
  var SampleSet = ee.FeatureCollection([TargetSample_cls,SourceSample_cls]).flatten();//,SourceSample_cls
  print(SourceSample_cls.size(),'SourceSample_cls');
  print(TargetSample_cls.size(),'TargetSample_cls');
  print(SampleSet.size(),'sampleset');
  
  var palette = ['#BC9F3D','#3A2700','#00755D','#00BCDB','#97B1AB','#F0DEB4']
  
  var classifier = ee.Classifier.smileRandomForest(130).train(SampleSet,'class',features.bandNames(),1,0); 
  print(classifier.explain());
  var classified = features.classify(classifier);
  Map.addLayer(classified,{palette:palette,max:6,min:1})
  
  var seeds = ee.Algorithms.Image.Segmentation.seedGrid(5);
  var snic = ee.Algorithms.Image.Segmentation.SNIC({
    image: features, 
    compactness: 0,
    connectivity: 4,
    neighborhoodSize: 10,
    size: 2,
    seeds: seeds
  })
  var clusters = snic.select('clusters')
  
  // Assign class to each cluster based on 'majority' voting (using ee.Reducer.mode()
  var smoothed = classified.addBands(clusters);
  
  var clusterMajority = smoothed.reduceConnectedComponents({
    reducer: ee.Reducer.mode(),
    labelBand: 'clusters'
  });
  Map.addLayer(clusterMajority, {palette:palette,max:6,min:1}, 
    'Processed using Clusters');
    
  //})
   Export.image.toDrive({
     image:clusterMajority,
     description:'cluster'+id,
     folder:'NCPClassified',
     fileNamePrefix :'cluster'+id,
     region:ROI,
     maxPixels:1e11,
     scale:10})
  Export.table.toDrive({
    collection:SourceSample11,
    description:'DD_TrSample'+id,
    folder:'TrSample',
    fileNamePrefix:'TrSample'+id});
       
}
