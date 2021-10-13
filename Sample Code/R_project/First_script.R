#This is my first code###
a=1+1
a=a+1
b="Lena"
b*a
b="3"

#######vector#####
a<-c(1,2,3)
b=c(2,3,4)
b=c(2,3,4,5)
##################
a+b
a-b
a*b
a%*%b
b%*%a
##################
a+b
a%*%b
##################
a[2]
b[1:2]
c=b[1:3]
d=b[-4]
###################
g=4
d=b[-g]
f=g*a
###################

v1<-c(a,b)
v1
v2<-c(1,2,3)
v2
v3=c(v1,4,5)
v3
v4=c(c(4,2,3),v3,c(3,8),1,2,c(v1,v2,v3))
v4
########Matrix#######


A=matrix(v4,2,16)

v4[3]
A[,1]
##################
Name=c("Anton","Berta","Caesar","Dora","Emil")
Gender=c("m","w","m","w","m")
Height=c(182,174,189,165,180)
Weight=c(80,68,92,55,78)
Name 
Gender
Height
Weight
##################
a=1:20
b=seq(1,20,le=20)
c=seq(1,20,by=1)
#####################
b=seq(1,20,le=100)
c=seq(1,20,by=0.1)
#######################
d=rep(1,10)
e=rep(a,3)
e=rep(a,each=3)
#######################
x<-(1:5)
a0<-sample(x,2,replace=2)
#set.seed(1)
a1<-sample(x,2)
a1
A <-rbind(a0,a1)
A
solve(A)
t(a1)
##################################

m1=matrix(c(1,4,3,0),2,2)
m2=matrix(c(1,2,3,2),2,2)
m1
m2
m1+1
m1+m2
m1-m2
m1*2
m1*m2
m1/2
m1/m2
m1%*%m2
t(m1)
solve(m1)
m1%*%solve(m1)

#Construction of lists
v=c("Hallo","Tag","Hi","Guten Tag")
m=matrix(1:9,3,3)
v
m
objects=list(v,m)
objects

#######################################
Name 
Gender
Height
Weight
#cbind(Name, Gender, Height, Weight)
Persons=data.frame(Name, Gender, Height, Weight)
summary(Persons)
Persons[1,2:3]
Persons$Name
Persons$Height

x=c(8,5,5,1,6,9,7,4)
y=1:length(x)
sum(x<=5)

