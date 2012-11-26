
#Their params:
A<-11.53
B<--1.78
C<-1833.599
D<-8834.987
E<-0.813
para<-c(A=A,B=B,C=C,D=D,E=E)

paraNames<-c('A','B','C','D','E')

#data:
df
dfty
dftxy #with 1/Y^2 & 1/Y^2

res<-drm(fi ~ expected_conc, data=df, fct=LL.5())
resty<-drm(fi ~ expected_conc, data=dfty, fct=LL.5())
restxy<-drm(fi ~ expected_conc, data=dftxy, fct=LL.5())

para<-res$coefficients
names(para)<-paraNames

func0<-function(x){para['A']+(para['D']/(1+(x/para['C'])^para['B'])^para['E'])}
plot(func0, from=1, log="xy", to=10000)