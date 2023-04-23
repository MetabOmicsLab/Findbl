##fooblet:xgboost+integration 
##内置人源和鼠源hashing数据，训练模型判断输入数据中s/d

##set parameters:输入数据以及参数(mod,hm或ms) and pred threshold
##输入数据格式

findbl <- function(count,eta=1,mod,pred=0.5){
  require(xgboost);require(caret)
  
  print("Training model...")
  
  # call inner model
  if (mod == "hm"){
    #xgbtrain after cv cv can set parameters and rounds
    hm <- readRDS("./train_data/human_integrated.rds")
    
    label <- hm$label
    data <- hm$counts
    n = nrow(data)
    train.index = sample(n,floor(0.75*n))
    
    train.data = data[train.index,]
    train.label = label[train.index]
    test.data = data[-train.index,]
    test.label = label[-train.index]
    
    dtrain <- xgb.DMatrix(data = train.data,label=train.label)
    dtest <- xgb.DMatrix(data = test.data,label=test.label)
    
    watchlist <- list(train=dtrain,test=dtest) 
    set.seed(123)
    
    bst <- xgb.train(data = dtrain,
                     nrounds = 3,
                     nthread = 20,
                     max_depth = 6,
                     eta = 1,
                     watchlist = watchlist,
                     eval.metric= "error",
                     eval.metric= "logloss",
                     eval.metric="auc",
                     objective = "binary:logistic")
    
    label<-getinfo(dtest,"label")
    pred<-predict(bst,dtest)
    xgb.save(bst,"hm_xgboost.model")
    boost_model <- xgb.load("./hm_xgboost.model")
  }
  
  if (mod == "ms"){
    ###
    ms <- readRDS("./train_data/mouse_multi_integrated.rds")
    
    label <- ms$label
    data <- ms$counts
    n = nrow(data)
    train.index = sample(n,floor(0.75*n))
    
    train.data = data[train.index,]
    train.label = label[train.index]
    test.data = data[-train.index,]
    test.label = label[-train.index]
    
    dtrain <- xgb.DMatrix(data = train.data,label=train.label)
    dtest <- xgb.DMatrix(data = test.data,label=test.label)
    
    watchlist <- list(train=dtrain,test=dtest) 
    set.seed(123)
    
    bst <- xgb.train(data = dtrain,
                     nrounds = 5,
                     nthread = 20,
                     max_depth = 6,
                     eta = 1,
                     watchlist = watchlist,
                     eval.metric= "error",
                     eval.metric= "logloss",
                     eval.metric="auc",
                     objective = "binary:logistic")
    
    label<-getinfo(dtest,"label")
    pred<-predict(bst,dtest)
    xgb.save(bst,"ms_xgboost.model")
    boost_model <- xgb.load("./ms_xgboost.model")
  }
  
  # use xgb model
  print("Predicting doublets...")

  ###count is seurat original count need transpose
  predict <- predict(boost_model,count)
  prediction <- as.numeric(predict>0.5)
  findbl.result <<- list(data = count, predict = predict,
                          prediction = prediction)
  print("Classification prediction complete.")
  
}

multi<-readRDS("./train_data/mouse_multi_integrated.rds")
findbl(count=multi$counts,mod="ms")

debug(findbl)
