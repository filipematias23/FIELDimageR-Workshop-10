############################
### FIELDimageR pipeline ###
############################

################
### Packages ### 
################

library(FIELDimageR)
library(raster)
library(ggplot2)
library(agricolae)
library(reshape2)
library(lme4)
library(readxl)

# Uploading file (odm_orthophoto.tif):
EX2<-stack("odm_orthophoto.tif")
plotRGB(EX2)

# Rotating the image with theta 2.3:
EX2.Rotated<-fieldRotate(mosaic = EX2[[-4]],theta = 2.3)

# Removing soil (example using HUE):
EX2.RemSoil<- fieldMask(mosaic = EX2.Rotated)

# Trial: 01
EX2.Crop.1 <- fieldCrop(mosaic = EX2.RemSoil$newMosaic) 

# Dataset - Field notes - phenotype (.exel or .txt):
Data1<-read_excel('EX2_Data.xlsx',1)

# Field map ID identification (ID="Plot"):
Map1<-fieldMap(fieldPlot = Data1$Plot,fieldColumn = Data1$Column, fieldRow = Data1$Row,decreasing = T)
Map1

# Building the plot shapefile (ncols = 14 and nrows = 10)
x11()
EX2.Shape.1<-fieldShape(mosaic = EX2.Crop.1, ncols = 14, nrows = 10, fieldData = Data1, ID = "Plot", fieldMap = Map1)

# Trial: 02
EX2.Crop.2 <- fieldCrop(mosaic = EX2.RemSoil$newMosaic) 
Data2<-read_excel('EX2_Data.xlsx',2)
Map2<-fieldMap(fieldPlot = Data2$Plot,fieldColumn = Data2$Column, fieldRow = Data2$Row,decreasing = T)
EX2.Shape.2<-fieldShape(mosaic = EX2.Crop.2, ncols = 14, nrows = 9, fieldData = Data2, ID = "Plot", fieldMap = Map2)

# Trial: 03
EX2.Crop.3 <- fieldCrop(mosaic = EX2.RemSoil$newMosaic) 
Data3<-read_excel('EX2_Data.xlsx',3)
Map3<-fieldMap(fieldPlot = Data3$Plot,fieldColumn = Data3$Column, fieldRow = Data3$Row,decreasing = T)
EX2.Shape.3<-fieldShape(mosaic = EX2.Crop.3, ncols = 14, nrows = 13, fieldData = Data3, ID = "Plot", fieldMap = Map3)

# Trial: 04
EX2.Crop.4 <- fieldCrop(mosaic = EX2.RemSoil$newMosaic) 
Data4<-read_excel('EX2_Data.xlsx',4)
Map4<-fieldMap(fieldPlot = Data4$Plot,fieldColumn = Data4$Column, fieldRow = Data4$Row,decreasing = T)
EX2.Shape.4<-fieldShape(mosaic = EX2.Crop.4, ncols = 14, nrows = 10, fieldData = Data4, ID = "Plot", fieldMap = Map4)

# Combining shapefiles
EX2.Shape<-rbind(EX2.Shape.1$fieldShape,EX2.Shape.2$fieldShape,EX2.Shape.3$fieldShape,EX2.Shape.4$fieldShape)
plot(EX2.Shape)

# Building indices ("NGRDI","BGI","GLI","SCI")
EX2.Indices<- fieldIndex(mosaic = EX2.RemSoil$newMosaic,
                         index = c("NGRDI","BGI"), 
                         myIndex = c("(Red-Blue)/Green"))

plot(EX2.Indices$myIndex)
plot(EX2.Shape,add=T)

# Extracting Data
EX2.Info<- fieldInfo(mosaic = EX2.Indices[[c("NGRDI","BGI","myIndex")]],
                     fieldShape = EX2.Shape,
                     buffer = -0.05,
                     n.core = 3)
NewData<-EX2.Info$fieldShape@data
NewData

# Making map plots
fieldPlot(fieldShape=EX2.Info$fieldShape,
          fieldAttribute="NGRDI", 
          mosaic=EX2.Indices, 
          color=c("red","blue"), 
          alpha = 0.4,
          round = 2)

# Uploading files from vegetative growth (dsm_1.tif):
DSM <- stack("dsm.tif")
plot(DSM)

# Rotating the image using the same theta (2.3):
DSM.R<-fieldRotate(DSM, theta = 2.3)

# Removing the soil using mask from step 4:
DSM.S <- fieldMask(DSM.R, mask = EX2.RemSoil$mask)

# Extracting the estimate plant height average (EPH):
EPH <- fieldInfo(DSM.S$newMosaic, fieldShape = EX2.Info$fieldShape, fun = "mean",n.core = 3)
EPH$plotValue
#EPH<-read.csv("EPH.csv",header = T) # In case of ERROR 

# Correlation
NewData2<-EPH$fieldShape@data[,c("dsm","NGRDI","BGI","myIndex")]
#NewData2<-EPH[,c("dsm","NGRDI","BGI","myIndex")]
cor1<-correlation(NewData2)
cor1$correlation
cor1$pvalue

# Regression
NewData3<-EPH$fieldShape@data[,c("Trial","Name","Block","Row","Column","dsm","NGRDI")]
#NewData3<-EPH[,c("Trial","Name","Block","Row","Column","dsm","NGRDI")]
NewData3$Trial<-as.factor(NewData3$Trial)
NewData3$Name<-as.factor(NewData3$Name)
NewData3$Block<-as.factor(NewData3$Block)
NewData3$Row<-as.factor(NewData3$Row)
NewData3$Column<-as.factor(NewData3$Column)
NewData3$dsm<-as.numeric(as.character(NewData3$dsm))
NewData3$NGRDI<-as.numeric(as.character(NewData3$NGRDI))

ggplot(NewData3,aes(x=dsm,y=NGRDI, col=Trial))+
  facet_wrap(~Trial, scales = "fixed",ncol = 1)+
  geom_point() +
  geom_smooth(method=lm)+
  theme_bw()

fieldPlot(fieldShape=EPH$fieldShape,
          fieldAttribute="dsm", 
          mosaic=EX2.Indices, 
          color=c("red","blue"), 
          alpha = 0.4,
          round = 2)

# Heritability - NGRDI

str(NewData3)

mod.T1<-lmer(NGRDI~Row+Column+(1|Name),NewData3[NewData3$Trial=="T1",])
H2.T1<-as.data.frame(VarCorr(mod.T1))$vcov[1]/sum(as.data.frame(VarCorr(mod.T1))$vcov)

mod.T2<-lmer(NGRDI~Row+Column+(1|Name),NewData3[NewData3$Trial=="T2",])
H2.T2<-as.data.frame(VarCorr(mod.T2))$vcov[1]/sum(as.data.frame(VarCorr(mod.T2))$vcov)

mod.T3<-lmer(NGRDI~Row+Column+(1|Name),NewData3[NewData3$Trial=="T3",])
H2.T3<-as.data.frame(VarCorr(mod.T3))$vcov[1]/sum(as.data.frame(VarCorr(mod.T3))$vcov)

mod.T4<-lmer(NGRDI~Row+Column+(1|Name),NewData3[NewData3$Trial=="T4",])
H2.T4<-as.data.frame(VarCorr(mod.T4))$vcov[1]/sum(as.data.frame(VarCorr(mod.T4))$vcov)

NewData4<-data.frame(Trail=factor(c("T1","T2","T3","T4")),
                     H2=as.numeric(c(H2.T1,H2.T2,H2.T3,H2.T4)))

ggplot(NewData4,aes(x=Trail,y=H2,fill=Trail))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(H2,2)), vjust=1.6, color="white", size=6)+
  theme_bw()

# Heritability - Estimated plant height 

mod.T1<-lmer(dsm~Row+Column+(1|Name),NewData3[NewData3$Trial=="T1",])
H2.T1<-as.data.frame(VarCorr(mod.T1))$vcov[1]/sum(as.data.frame(VarCorr(mod.T1))$vcov)

mod.T2<-lmer(dsm~Row+Column+(1|Name),NewData3[NewData3$Trial=="T2",])
H2.T2<-as.data.frame(VarCorr(mod.T2))$vcov[1]/sum(as.data.frame(VarCorr(mod.T2))$vcov)

mod.T3<-lmer(dsm~Row+Column+(1|Name),NewData3[NewData3$Trial=="T3",])
H2.T3<-as.data.frame(VarCorr(mod.T3))$vcov[1]/sum(as.data.frame(VarCorr(mod.T3))$vcov)

mod.T4<-lmer(dsm~Row+Column+(1|Name),NewData3[NewData3$Trial=="T4",])
H2.T4<-as.data.frame(VarCorr(mod.T4))$vcov[1]/sum(as.data.frame(VarCorr(mod.T4))$vcov)

NewData4<-data.frame(Trail=factor(c("T1","T2","T3","T4")),
                     EPH=as.numeric(c(H2.T1,H2.T2,H2.T3,H2.T4)))

ggplot(NewData4,aes(x=Trail,y=EPH,fill=Trail))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(EPH,2)), vjust=1.6, color="white", size=6)+
  theme_bw()

###########
### END ###
###########
