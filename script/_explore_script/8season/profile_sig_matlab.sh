#! /bin/sh
octave <<mark
pkg load signal
close all
clear all
control=load ('-ascii','$1');
sensi=load ('-ascii','$2');
alpha=0.05;
col1=control(:,1);
col2=control(:,2);
col3=sensi(:,1);
col4=sensi(:,2);
n=length(col2);
for i=1:10000
sampy1=col2(int32(rand(n,1)*(n-1))+1);
%sampy1=col2(ceil(rand(n,1)*n);
sampy2=col4(int32(rand(n,1)*(n-1))+1);
bootstrap(i)=mean(sampy1)-mean(sampy2);
end
confidence=prctile(bootstrap,[100*alpha/2,100-alpha/2]);
significance=confidence(1)>0 | confidence(2) <0;
save ('-ascii','${1}_${2}_significant_different_mean_5_percent_confidence_level','confidence','significance');
disp('confidence interval');
disp(confidence);
disp('Does the confidence interval cover zero (0=yes; 1=no)');
disp(significance);
mark
