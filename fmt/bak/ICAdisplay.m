./ViewFMRI.m                                                                                        000644  002226  001747  00000020133 07433525376 014446  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [tseries,outarg] = ViewFMRI(Data,Mean,func,funcd,varargin)

% tseries = ViewFMRI(Data)
%
% Tool for viewing 4D FMRI datasets
%
% tseries = ViewFMRI(Data,Mean)
%
% Mean determines the data shown in the slice-by-slice lightbox.
% Mean='m' shows the mean volume, this is the default.
% Can alternatively be a matrix volume of the same size as the data
%
% tseries = ViewFMRI(Data,Mean,func,funcd,varargin)
%
% func is function that runs on middle mouse click
% default is plotseries, which just plots the time series
% an alternative is fftseries, which plots the fft of the time series
%
% funcd is function that runs when d is pressed
% default is plotseries, which just plots the time series
% an alternative is fftseries, which plots the fft of the time series
%
% varargin are any arguments to be passed to func and funcd

if (nargin<4)
  funcd = 'plotseries';
  if (nargin<3)
    func = 'plotseries';
  end;
end;

outarg=0;

if size(Data,1) ~= size(Data,2)
  fprintf('Image matrix has been padded with zeros to make it square');
  d1 = max([size(Data,1),size(Data,2)]);
  Data2 = zeros(d1,d1,size(Data,3),size(Data,4));
  Data2(1:size(Data,1),1:size(Data,2),:,:) = Data;
  Data = Data2;
  clear Data2;
end;

res = size(Data,1);

Index=0;

if (nargin<2),
  Mean='m';
end;

if length(size(Data))==2
  NumSlices=1;
end;

% Beckmanns old stuff for viewing 2D matrix
% if length(size(Data))==2
%  MeanImg=Data(floor(size(Data,1)/2),:);
%  if (nargin>1)
%    if Mean=='m'
%      MeanImg=mean(Data);
%    end;
%  end;
%  
%  NumSlices=size(Data,2)/slicesize;
%  Data=reshape(Data',res,res,NumSlices,size(Data,1));
%end;


if ~isstr(Mean),
  MeanImg=Mean;

  if size(MeanImg,1) ~= size(MeanImg,2)
      d1 = max([size(MeanImg,1),size(MeanImg,2)]);
      MeanImg2 = zeros(d1,d1,size(MeanImg,3));
      MeanImg2(1:size(MeanImg,1),1:size(MeanImg,2),:) = MeanImg;
      MeanImg = MeanImg2;
      clear MeanImg2;
  end;

else,
  MeanImg=mean(Data,4);
end;


NumSlices=size(MeanImg,3);

if NumSlices>21
      AllData = MeanImg;
      MeanImg = AllData(:,:,Index+1:Index+21); 
end
MeanImg=reshape(MeanImg,1,prod(size(MeanImg)));

if NumSlices<21
  tmp=zeros(1,size(Data,1)*size(Data,2)*21);
  tmp(1:size(Data,1)*size(Data,2)*NumSlices)=MeanImg;
  MeanImg=tmp;
end;

MeanImg=reshape(MeanImg,res,res,21);

LeftPanel=zeros(7*res,3*res);
for ct1=0:2
  for ct2=0:6
    LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)=flipud(MeanImg(:,:,1+ct1+ct2*3)');
  end;
end;

HndlImg=figure;

factor=2;

colormap('bone');
Xrange=[1,res];
Yrange=[1,res];
scale=res;
ImgNum=1;Xind=1;Yind=1;  

LeftImg=subplot('position',[0.01,0.05,0.29,0.9]);
imagesc(LeftPanel);
axis('equal')
axis('off')
if NumSlices > 21
  title(sprintf('n for next 21 slices\np for previous 21 slices'));
end;

tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
MeanInt=MeanImg(Xind,Yind,ImgNum+Index);
subplot('position',[0.37,0.05,0.61,0.29]);    
feval(func, tseries, Xind,Yind,ImgNum+Index,varargin{:});

RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
imagesc((MeanImg(:,:,ImgNum)'))
set(RightImg,'XLim',Xrange,'YLim',Yrange);
axis('equal')
axis('tight')
set(RightImg,'XColor','red','YColor','red');


while 1==1
  
  TextImg=subplot('position',[0.8,0.41,0.18,0.55]);
  plot([1])
  axis('off')
  set(TextImg,'XLim',[0,2],'YLim',[0,2]);
  text(0.3,1.8,sprintf('Slice No: %2d',ImgNum+Index));
  text(0.3,1.6,sprintf('       X: %2d',Xind));
  text(0.3,1.4,sprintf('       Y: %2d',Yind));
  text(0.3,1.2,sprintf('Intensity'));
  text(0.3,1.0,sprintf('     %5.3f',MeanInt));
  line([0.05,1.95],[0.85,0.85],'Color','black')
  text(0.05,0.7,sprintf('<left>,<right>: zoom'));
  text(0.05,0.5,sprintf('<centre>: plot tc '));
  text(0.05,0.3,sprintf('q,Q: quit'));
  text(0.05,0.1,sprintf('d,D: indiv. tc'));
  
  drawflag=1;
  
  subplot(RightImg); 
  [X_c,Y_c,B_c]=ginput(1);
  
  if gca==LeftImg	 
    
    tmp=1+floor(X_c/res)+3*floor(Y_c/res);
    if tmp<=NumSlices-Index
      ImgNum=tmp;
    else
      drawflag=0;
    end;
    
    if drawflag
      tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
      
      MeanInt=MeanImg(Xind,Yind,ImgNum);
      subplot('position',[0.37,0.05,0.61,0.29]);    
      feval(func, tseries, Xind,Yind,ImgNum+Index,varargin{:});
      
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc((MeanImg(:,:,ImgNum)'))
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum+Index)));
      set(RightImg,'XColor','red','YColor','red');
      
      LeftImg=subplot('position',[0.01,0.05,0.29,0.9]);
      imagesc(LeftPanel);
      axis('off')
      axis('equal')
      Corner1=1+res*floor(X_c/res);
      Corner2=1+res*floor(Y_c/res);
      line([Corner1,Corner1],[Corner2,Corner2+res-1]);
      line([Corner1+res-1,Corner1+res-1],[Corner2,Corner2+res-1]);
      line([Corner1,Corner1+res-1],[Corner2,Corner2]);
      line([Corner1,Corner1+res-1],[Corner2+res-1,Corner2+res-1]);
    end;
  end; 
  
  if gca==RightImg
    if (B_c==1)|(B_c==2)|(B_c==3) 
      if (round(X_c)>0)&(round(X_c)<=res)
	Xind=round(X_c);
      end;
      if (round(Y_c)>0)&(round(Y_c)<=res)      
	Yind=round(Y_c);
      end;
    end;
    if B_c==2
      TSImg=subplot('position',[0.37,0.05,0.61,0.29]);    
      tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
      MeanInt=MeanImg(Xind,Yind,ImgNum);
      feval(func, tseries, Xind,Yind,ImgNum+Index,varargin{:});
    end;
    if (B_c==1)&(scale>3)
      scale=round(scale/factor);
      Xrange=[Xind-0.5*scale,Xind+0.5*scale-1];
      Yrange=[Yind-0.5*scale,Yind+0.5*scale-1];
      %%Xrange=Xrange+sum(Xrange<1)-sum(Xrange>res);
      %%Yrange=Yrange+sum(Yrange<1)-sum(Yrange>res);
      Xrange=Xrange-min(Xrange-1)*max(Xrange<1)-max(Xrange-res)*max(Xrange>res);
      Yrange=Yrange-min(Yrange-1)*max(Yrange<1)-max(Yrange-res)*max(Yrange>res);
      %%fprintf('X:%d Y:%d Xr: %d-%d Yr: %d-%d\n',Xind,Yind,min(Xrange),max(Xrange),min(Yrange),max(Yrange))
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc((MeanImg(:,:,ImgNum)'))
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum+Index)));
      set(RightImg,'XColor','red','YColor','red');
    end;
    if (B_c==3)&(scale < res-1)
      scale=round(scale*factor);
      Xrange=[Xind-0.5*scale,Xind+0.5*scale-1];
      Yrange=[Yind-0.5*scale,Yind+0.5*scale-1];
      Xrange=Xrange-min(Xrange-1)*max(Xrange<1)-max(Xrange-res)*max(Xrange>res);
      Yrange=Yrange-min(Yrange-1)*max(Yrange<1)-max(Yrange-res)*max(Yrange>res);
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc((MeanImg(:,:,ImgNum)'))
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum+Index)));
      set(RightImg,'XColor','red','YColor','red');
    end;
  end;
  
  if (B_c==113)|(B_c==81)
    % close(HndlImg);
    break
  end;	 
  if (B_c==68)|(B_c==100)
    figure
    feval(funcd, tseries, Xind,Yind,ImgNum+Index,varargin{:});     
    if size(outarg,1)==1
       outarg=tseries;
    else
       outarg=[outarg,tseries];
    end;
    
    
    figure(HndlImg);
  end;
  if (B_c==112)|(B_c==80)
    %P or p
    if (Index>0)
      Index=Index-21;
      MeanImg = zeros(res,res,21);
      MeanImg(:,:,1:min([Index+21 NumSlices])-Index) = AllData(:,:,Index+1:min([Index+21 NumSlices])); 
      LeftPanel=zeros(7*res,3*res);
      for ct1=0:2
	for ct2=0:6
	  LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)=MeanImg(:,:,1+ct1+ct2*3)';
	end;
      end;
      subplot(LeftImg);
      imagesc(LeftPanel);
      axis('equal')
      axis('off')
      if NumSlices > 21
	title(sprintf('n for next 21 slices\np for previous 21 slices'));
      end;
    end;
  end;
  if (B_c==110)|(B_c==78)
    %N or n
    if ((Index+21)<NumSlices)
      Index=Index+21;
      MeanImg = zeros(res,res,21);
      MeanImg(:,:,1:min([Index+21 NumSlices])-Index) = AllData(:,:,Index+1:min([Index+21 NumSlices])); 
      LeftPanel=zeros(7*res,3*res);
      for ct1=0:2
	for ct2=0:6
	  LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)=MeanImg(:,:,1+ct1+ct2*3)';
	end;
      end;
      subplot(LeftImg);
      imagesc(LeftPanel);
      axis('equal')
      axis('off')
      if NumSlices > 21
	title(sprintf('n for next 21 slices\np for previous 21 slices'));
      end;
    end;
  end;
end;
                                                                                                                                                                                                                                                                                                                                                                                                                                     ./ViewFMRI2.m                                                                                       000664  002226  001747  00000017631 07323360060 014524  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function tseries = ViewFMRI2(Data,Mean,func,funcd,xdim,ydim,varargin)

% tseries = ViewFMRI(Data,Mean,func,funcd,varargin)
%
% Tool for viewing 4D FMRI datasets.
%
% Mean determines the data shown in the slice-by-slice lightbox.
% Mean='m' shows the mean volume.
% Mean='' shows the middle volume in the series, this is the
% default.
%
% func is function that runs on middle mouse click
% default is plotseries, which just plots the time series
% an alternative is fftseries, which plots the fft of the time series
%
% funcd is function that runs when d is pressed
% default is plotseries, which just plots the time series
% an alternative is fftseries, which plots the fft of the time series
%
% varargin are any arguments to be passed to func and funcd
%
% xdim,ydim are voxel dimensions - default 64,64

   
if (nargin<6)
   ydim=64;
   if (nargin<5)
      xdim=64;
      if (nargin<4)
         funcd = 'plotseries';
         if (nargin<3)
            func = 'plotseries';
	 end;
      end;
   end;
end;

%if size(Data,1) ~= size(Data,2)
%  error('Image matrix must be square');
%end;
%res = size(Data,1);

Index=0;
if length(size(Data))==2
  MeanImg=Data(floor(size(Data,1)/2),:);
  if (nargin>1)
    if Mean=='m'
      MeanImg=mean(Data);
    end;
  end;
  
  NumSlices=size(Data,2)/(xdim*ydim);
  Data=reshape(Data',xdim,ydim,NumSlices,size(Data,1));

else
  MeanImg=mean(Data,4);
  NumSlices=size(MeanImg,3);
  if NumSlices>21
    AllData = MeanImg;
    MeanImg = AllData(:,:,Index+1:Index+21); 
  end
  MeanImg=reshape(MeanImg,1,prod(size(MeanImg)));
end;

if NumSlices<21
  tmp=zeros(1,xdim*ydim*21);
  tmp(1:xdim*ydim*NumSlices)=MeanImg;
  MeanImg=tmp;
end;

MeanImg=reshape(MeanImg,xdim,ydim,21);

LeftPanel=zeros(7*ydim,3*xdim);
for ct1=0:2
  for ct2=0:6
    LeftPanel(1+ct2*ydim:ydim+ct2*ydim,1+ct1*xdim:xdim+ct1*xdim)=MeanImg(:,:,1+ct1+ct2*3)';
  end;
end;


HndlImg=figure;

factor=2;

colormap('bone');
Xrange=[1,xdim];
Yrange=[1,ydim];
scale=min(xdim,ydim);
ImgNum=1;Xind=1;Yind=1;  

LeftImg=subplot('position',[0.01,0.05,0.29,0.9]);
imagesc(LeftPanel);
if xdim==ydim
   axis('equal')
end;
axis('off')
if NumSlices > 21
  title(sprintf('n for next 21 slices\np for previous 21 slices'));
end;

tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
MeanInt=mean(tseries);
subplot('position',[0.37,0.05,0.61,0.29]);    
feval(func, tseries, varargin{:});
title(sprintf('Slice %d; X: %d  Y: %d -- Mean Intensity: %5.3f ',ImgNum-1+Index,Xind,Yind,MeanInt)); 

RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
imagesc(MeanImg(:,:,ImgNum)')
set(RightImg,'XLim',Xrange,'YLim',Yrange);
if xdim==ydim
   axis('equal')
end;
axis('tight')
set(RightImg,'XColor','red','YColor','red');


while 1==1
  
  TextImg=subplot('position',[0.8,0.41,0.18,0.55]);
  plot([1])
  axis('off')
  set(TextImg,'XLim',[0,2],'YLim',[0,2]);
  text(0.3,1.8,sprintf('Slice No: %2d',ImgNum-1+Index));
  text(0.3,1.6,sprintf('       X: %2d',Xind));
  text(0.3,1.4,sprintf('       Y: %2d',Yind));
  text(0.3,1.2,sprintf('Mean Intensity'));
  text(0.3,1.0,sprintf('     %5.3f',MeanInt));
  line([0.05,1.95],[0.85,0.85],'Color','black')
  text(0.05,0.7,sprintf('<left>,<right>: zoom'));
  text(0.05,0.5,sprintf('<centre>: plot tc '));
  text(0.05,0.3,sprintf('q,Q: quit'));
  text(0.05,0.1,sprintf('d,D: indiv. tc'));
  
  drawflag=1;
  
  subplot(RightImg); 
  [X_c,Y_c,B_c]=ginput(1);
  
  if gca==LeftImg	 
    
    tmp=1+floor(X_c/xdim)+3*floor(Y_c/ydim);
    if tmp<=NumSlices-Index
      ImgNum=tmp;
    else
      drawflag=0;
    end;
    
    if drawflag
      tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
      MeanInt=mean(tseries);
      subplot('position',[0.37,0.05,0.61,0.29]);    
      feval(func, tseries, varargin{:});
      title(sprintf('Slice %d; X: %d  Y: %d -- Mean Intensity: %5.3f ',ImgNum-1+Index,Xind,Yind,MeanInt)); 
      
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc(MeanImg(:,:,ImgNum)')
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum-1+Index)));
      set(RightImg,'XColor','red','YColor','red');
      
      LeftImg=subplot('position',[0.01,0.05,0.29,0.9]);
      imagesc(LeftPanel);
      axis('off')
      if xdim==ydim
	 axis('equal')
      end;
      
      Corner1=1+xdim*floor(X_c/xdim);
      Corner2=1+ydim*floor(Y_c/ydim);
      line([Corner1,Corner1],[Corner2,Corner2+ydim-1]);
      line([Corner1+xdim-1,Corner1+xdim-1],[Corner2,Corner2+ydim-1]);
      line([Corner1,Corner1+xdim-1],[Corner2,Corner2]);
      line([Corner1,Corner1+xdim-1],[Corner2+ydim-1,Corner2+ydim-1]);
    end;
  end; 
  
  if gca==RightImg
    if (B_c==1)|(B_c==2)|(B_c==3) 
      if (round(X_c)>0)&(round(X_c)<=xdim)
	Xind=round(X_c);
      end;
      if (round(Y_c)>0)&(round(Y_c)<=ydim)      
	Yind=round(Y_c);
      end;
    end;
    if B_c==2
      TSImg=subplot('position',[0.37,0.05,0.61,0.29]);    
      tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
      MeanInt=mean(tseries);
      feval(func, tseries, varargin{:});
      %title(sprintf('Slice %d; X: %d  Y: %d -- Mean Intensity: %5.3f ',ImgNum-1+Index,Xind,Yind,MeanInt)); 
    end;
    if (B_c==1)&(scale>3)
      scale=round(scale/factor);
      Xrange=[Xind-0.5*scale,Xind+0.5*scale-1];
      Yrange=[Yind-0.5*scale,Yind+0.5*scale-1];
      %%Xrange=Xrange+sum(Xrange<1)-sum(Xrange>xdim);
      %%Yrange=Yrange+sum(Yrange<1)-sum(Yrange>ydim);
      Xrange=Xrange-min(Xrange-1)*max(Xrange<1)-max(Xrange-xdim)*max(Xrange>xdim);
      Yrange=Yrange-min(Yrange-1)*max(Yrange<1)-max(Yrange-ydim)*max(Yrange>ydim);
      %%fprintf('X:%d Y:%d Xr: %d-%d Yr: %d-%d\n',Xind,Yind,min(Xrange),max(Xrange),min(Yrange),max(Yrange))
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc(MeanImg(:,:,ImgNum)')
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum-1+Index)));
      set(RightImg,'XColor','red','YColor','red');
    end;
    if (B_c==3)&(scale < min(xdim,ydim)-1)
      scale=round(scale*factor);
      Xrange=[Xind-0.5*scale,Xind+0.5*scale-1];
      Yrange=[Yind-0.5*scale,Yind+0.5*scale-1];
      Xrange=Xrange-min(Xrange-1)*max(Xrange<1)-max(Xrange-xdim)*max(Xrange>xdim);
      Yrange=Yrange-min(Yrange-1)*max(Yrange<1)-max(Yrange-ydim)*max(Yrange>ydim);
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc(MeanImg(:,:,ImgNum)');
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum-1+Index)));
      set(RightImg,'XColor','red','YColor','red');
    end;
  end;
  
  if (B_c==113)|(B_c==81)
    close(HndlImg);
    break
  end;	 
  if (B_c==68)|(B_c==100)
    figure
    feval(funcd, tseries, varargin{:});     
    %title(sprintf('Slice %d; X: %d  Y: %d -- Mean Intensity: %5.3f ',ImgNum-1+Index,Xind,Yind,MeanInt));
    figure(HndlImg);
  end;
  if (B_c==112)|(B_c==80)
    %P or p
    if (Index>0)
      Index=Index-21;
      MeanImg = zeros(xdim,ydim,21);
      MeanImg(:,:,1:min([Index+21 NumSlices])-Index) = AllData(:,:,Index+1:min([Index+21 NumSlices])); 
      LeftPanel=zeros(7*ydim,3*xdim);
      for ct1=0:2
	for ct2=0:6
	  LeftPanel(1+ct2*ydim:ydim+ct2*ydim,1+ct1*xdim:xdim+ct1*xdim)=MeanImg(:,:,1+ct1+ct2*3)';
	end;
      end;
      subplot(LeftImg);
      imagesc(LeftPanel);
      if xdim==ydim
         axis('equal')
      end;
      axis('off')
      if NumSlices > 21
	title(sprintf('n for next 21 slices\np for previous 21 slices'));
      end;
    end;
  end;
  if (B_c==110)|(B_c==78)
    %N or n
    if ((Index+21)<NumSlices)
      Index=Index+21;
      MeanImg = zeros(xdim,ydim,21);
      MeanImg(:,:,1:min([Index+21 NumSlices])-Index) = AllData(:,:,Index+1:min([Index+21 NumSlices])); 
      LeftPanel=zeros(7*ydim,3*xdim);
      for ct1=0:2
	for ct2=0:6
	  LeftPanel(1+ct2*xdim:xdim+ct2*xdim,1+ct1*ydim:ydim+ct1*ydim)=MeanImg(:,:,1+ct1+ct2*3)';
	end;
      end;
      subplot(LeftImg);
      imagesc(LeftPanel);
      if xdim==ydim
         axis('equal')
      end;
      axis('off')
      if NumSlices > 21
	title(sprintf('n for next 21 slices\np for previous 21 slices'));
      end;
    end;
  end;
end;
                                                                                                       ./ViewFMRI3.m                                                                                       000644  002226  001747  00000020025 07351615365 014526  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function tseries = ViewFMRI3(Data,Mean,func,funcd,varargin)

% tseries = ViewFMRI(Data)
%
% Tool for viewing 4D FMRI datasets
%
% tseries = ViewFMRI(Data,Mean)
%
% Mean determines the data shown in the slice-by-slice lightbox.
% Mean='m' shows the mean volume, this is the default.
% Can alternatively be a matrix volume of the same size as the data
%
% tseries = ViewFMRI(Data,Mean,func,funcd,varargin)
%
% func is function that runs on middle mouse click
% default is plotseries, which just plots the time series
% an alternative is fftseries, which plots the fft of the time series
%
% funcd is function that runs when d is pressed
% default is plotseries, which just plots the time series
% an alternative is fftseries, which plots the fft of the time series
%
% varargin are any arguments to be passed to func and funcd

if (nargin<4)
  funcd = 'plotseries';
  if (nargin<3)
    func = 'plotseries';
  end;
end;

if size(Data,1) ~= size(Data,2)
  fprintf('Image matrix has been padded with zeros to make it square');
  d1 = max([size(Data,1),size(Data,2)]);
  Data2 = zeros(d1,d1,size(Data,3),size(Data,4));
  Data2(1:size(Data,1),1:size(Data,2),:,:) = Data;
  Data = Data2;
  clear Data2;
end;

res = size(Data,1);

Index=0;

if (nargin<2),
  Mean='m';
end;

if length(size(Data))==2
  NumSlices=1;
end;

% Beckmanns old stuff for viewing 2D matrix
% if length(size(Data))==2
%  MeanImg=Data(floor(size(Data,1)/2),:);
%  if (nargin>1)
%    if Mean=='m'
%      MeanImg=mean(Data);
%    end;
%  end;
%  
%  NumSlices=size(Data,2)/slicesize;
%  Data=reshape(Data',res,res,NumSlices,size(Data,1));
%end;


if ~isstr(Mean),
  MeanImg=Mean;

  if size(MeanImg,1) ~= size(MeanImg,2)
      d1 = max([size(MeanImg,1),size(MeanImg,2)]);
      MeanImg2 = zeros(d1,d1,size(MeanImg,3));
      MeanImg2(1:size(MeanImg,1),1:size(MeanImg,2),:) = MeanImg;
      MeanImg = MeanImg2;
      clear MeanImg2;
  end;

else,
  MeanImg=mean(Data,4);
end;


NumSlices=size(MeanImg,3);

if NumSlices>21
      AllData = MeanImg;
      MeanImg = AllData(:,:,Index+1:Index+21); 
end
MeanImg=reshape(MeanImg,1,prod(size(MeanImg)));

if NumSlices<21
  tmp=zeros(1,size(Data,1)*size(Data,2)*21);
  tmp(1:size(Data,1)*size(Data,2)*NumSlices)=MeanImg;
  MeanImg=tmp;
end;

MeanImg=reshape(MeanImg,res,res,21);

LeftPanel=zeros(7*res,3*res);
for ct1=0:2
  for ct2=0:6
    LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)=flipud(MeanImg(:,:,1+ct1+ct2*3)');
  end;
end;

HndlImg=figure;

factor=2;

colormap('bone');
Xrange=[1,res];
Yrange=[1,res];
scale=res;
ImgNum=1;Xind=1;Yind=1;  

LeftImg=subplot('position',[0.01,0.05,0.29,0.9]);
imagesc(LeftPanel);
axis('equal')
axis('off')
if NumSlices > 21
  title(sprintf('n for next 21 slices\np for previous 21 slices'));
end;

tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
MeanInt=MeanImg(Xind,Yind,ImgNum+Index);
subplot('position',[0.37,0.05,0.61,0.29]);    
feval(func, tseries, Xind,Yind,ImgNum+Index,varargin{:});

RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
imagesc((MeanImg(:,:,ImgNum)'))
set(RightImg,'XLim',Xrange,'YLim',Yrange);
axis('equal')
axis('tight')
set(RightImg,'XColor','red','YColor','red');
axis xy

while 1==1
  
  TextImg=subplot('position',[0.8,0.41,0.18,0.55]);
  plot([1])
  axis('off')
  set(TextImg,'XLim',[0,2],'YLim',[0,2]);
  text(0.3,1.8,sprintf('Slice No: %2d',ImgNum+Index));
  text(0.3,1.6,sprintf('       X: %2d',Xind));
  text(0.3,1.4,sprintf('       Y: %2d',Yind));
  text(0.3,1.2,sprintf('Intensity'));
  text(0.3,1.0,sprintf('     %5.3f',MeanInt));
  line([0.05,1.95],[0.85,0.85],'Color','black')
  text(0.05,0.7,sprintf('<left>,<right>: zoom'));
  text(0.05,0.5,sprintf('<centre>: plot tc '));
  text(0.05,0.3,sprintf('q,Q: quit'));
  text(0.05,0.1,sprintf('d,D: indiv. tc'));
  
  drawflag=1;
  
  subplot(RightImg); 
  [X_c,Y_c,B_c]=ginput(1);
  
  if gca==LeftImg	 
    
    tmp=1+floor(X_c/res)+3*floor(Y_c/res);
    if tmp<=NumSlices-Index
      ImgNum=tmp;
    else
      drawflag=0;
    end;
    
    if drawflag
      tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
      
      MeanInt=MeanImg(Xind,Yind,ImgNum);
      subplot('position',[0.37,0.05,0.61,0.29]);    
      feval(func, tseries, Xind,Yind,ImgNum+Index,varargin{:});
      
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc((MeanImg(:,:,ImgNum)'))
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum+Index)));
      set(RightImg,'XColor','red','YColor','red');
      axis xy;
      
      LeftImg=subplot('position',[0.01,0.05,0.29,0.9]);
      imagesc(LeftPanel);
      axis('off')
      axis('equal')
      Corner1=1+res*floor(X_c/res);
      Corner2=1+res*floor(Y_c/res);
      line([Corner1,Corner1],[Corner2,Corner2+res-1]);
      line([Corner1+res-1,Corner1+res-1],[Corner2,Corner2+res-1]);
      line([Corner1,Corner1+res-1],[Corner2,Corner2]);
      line([Corner1,Corner1+res-1],[Corner2+res-1,Corner2+res-1]);
    end;
  end; 
  
  if gca==RightImg
    if (B_c==1)|(B_c==2)|(B_c==3) 
      if (round(X_c)>0)&(round(X_c)<=res)
	Xind=round(X_c);
      end;
      if (round(Y_c)>0)&(round(Y_c)<=res)      
	Yind=round(Y_c);
      end;
    end;
    if B_c==2
      TSImg=subplot('position',[0.37,0.05,0.61,0.29]);    
      tseries=squeeze(Data(Xind,Yind,ImgNum+Index,:));
      MeanInt=MeanImg(Xind,Yind,ImgNum);
      feval(func, tseries, Xind,Yind,ImgNum+Index,varargin{:});
    end;
    if (B_c==1)&(scale>3)
      scale=round(scale/factor);
      Xrange=[Xind-0.5*scale,Xind+0.5*scale-1];
      Yrange=[Yind-0.5*scale,Yind+0.5*scale-1];
      %%Xrange=Xrange+sum(Xrange<1)-sum(Xrange>res);
      %%Yrange=Yrange+sum(Yrange<1)-sum(Yrange>res);
      Xrange=Xrange-min(Xrange-1)*max(Xrange<1)-max(Xrange-res)*max(Xrange>res);
      Yrange=Yrange-min(Yrange-1)*max(Yrange<1)-max(Yrange-res)*max(Yrange>res);
      %%fprintf('X:%d Y:%d Xr: %d-%d Yr: %d-%d\n',Xind,Yind,min(Xrange),max(Xrange),min(Yrange),max(Yrange))
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc((MeanImg(:,:,ImgNum)'))
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum+Index)));
      set(RightImg,'XColor','red','YColor','red');
      axis xy;
    end;
    if (B_c==3)&(scale < res-1)
      scale=round(scale*factor);
      Xrange=[Xind-0.5*scale,Xind+0.5*scale-1];
      Yrange=[Yind-0.5*scale,Yind+0.5*scale-1];
      Xrange=Xrange-min(Xrange-1)*max(Xrange<1)-max(Xrange-res)*max(Xrange>res);
      Yrange=Yrange-min(Yrange-1)*max(Yrange<1)-max(Yrange-res)*max(Yrange>res);
      RightImg=subplot('position',[0.37,0.41,0.4,0.55]);
      imagesc(flipud(MeanImg(:,:,ImgNum)'))
      set(RightImg,'XLim',Xrange,'YLim',Yrange);
      title(strcat('Slice ',num2str(ImgNum+Index)));
      set(RightImg,'XColor','red','YColor','red');
      axis xy
    end;
  end;
  
  if (B_c==113)|(B_c==81)
    % close(HndlImg);
    break
  end;	 
  if (B_c==68)|(B_c==100)
    figure
    feval(funcd, tseries, Xind,Yind,ImgNum+Index,varargin{:});     
    figure(HndlImg);
  end;
  if (B_c==112)|(B_c==80)
    %P or p
    if (Index>0)
      Index=Index-21;
      MeanImg = zeros(res,res,21);
      MeanImg(:,:,1:min([Index+21 NumSlices])-Index) = AllData(:,:,Index+1:min([Index+21 NumSlices])); 
      LeftPanel=zeros(7*res,3*res);
      for ct1=0:2
	for ct2=0:6
	  LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)=MeanImg(:,:,1+ct1+ct2*3)';
	end;
      end;
      subplot(LeftImg);
      imagesc(LeftPanel);
      axis('equal')
      axis('off')
      if NumSlices > 21
	title(sprintf('n for next 21 slices\np for previous 21 slices'));
      end;
    end;
  end;
  if (B_c==110)|(B_c==78)
    %N or n
    if ((Index+21)<NumSlices)
      Index=Index+21;
      MeanImg = zeros(res,res,21);
      MeanImg(:,:,1:min([Index+21 NumSlices])-Index) = AllData(:,:,Index+1:min([Index+21 NumSlices])); 
      LeftPanel=zeros(7*res,3*res);
      for ct1=0:2
	for ct2=0:6
	  LeftPanel(1+ct2*res:res+ct2*res,1+ct1*res:res+ct1*res)=MeanImg(:,:,1+ct1+ct2*3)';
	end;
      end;
      subplot(LeftImg);
      imagesc(LeftPanel);
      axis('equal')
      axis('off')
      if NumSlices > 21
	title(sprintf('n for next 21 slices\np for previous 21 slices'));
      end;
    end;
  end;
end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ./acseries.m                                                                                        000644  002206  001747  00000000273 07434762404 014721  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function acseries(tseries, i, j, k, varargin)

tseries = detrend(tseries,0);
ac = xcorr(tseries,20,'coeff');

plot(ac(ceil(end/2):end));

title(sprintf('Slice %d; X: %d  Y: %d',k,i,j)); 
                                                                                                                                                                                                                                                                                                                                     ./applyaffine.m                                                                                     000664  002226  001747  00000001344 07323360060 015342  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [warpedim] = applyaffine(im,aff)
% [warpedim] = applyaffine(im,aff)
%
% Applies the affine transform (aff) to image im, producing a 
%  transformed image (warpedim)
% The affine matrix must be either 3x3 or 4x4 but always using
%   homogeneous coordinates  (the 4x4 describes a 3D transform)

warpedim = zeros(size(im));
[N,M] = size(im);

if size(aff) == [3,3],
 aff3 = aff;
elseif size(aff) == [4,4],
 aff3 = [ aff(1:2,1:2) aff(1:2,4) ; 0 0 1];
else
 error('Affine matrix is not 3x3 or 4x4');
end

x=kron(1:N,ones(M,1));
y=kron((1:M)',ones(1,N));
aff3=aff3-eye(3);

warp = zeros(N,M,2);
warp(:,:,1) = aff3(1,1)*x + aff3(1,2)*y + aff3(1,3);
warp(:,:,2) = aff3(2,1)*x + aff3(2,2)*y + aff3(2,3);

warpedim = applywarp(im,warp,x,y);
                                                                                                                                                                                                                                                                                            ./applywarp.m                                                                                       000664  002226  001747  00000001327 07323360060 015064  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [warpedim] = applywarp(im,warp,x,y)
% [warpedim] = applywarp(im,warp,x,y)
%
% Applies the warping field (warp) to image im, producing a new
%  warped image
% The warp field must be size [N,M,2] where warp(:,:,1) is the
%  warp in x, and warp(:,:,2) is the warp in y
% The additional arguments, x and y, must give the z and y coordinates
%  (respectively) at each pixel : they need to be the same size as im

warpedim = zeros(size(im));
[N,M] = size(im);

wx = round(x + warp(:,:,1)); 
wx=min(max(wx,1),N);
wy = round(y + warp(:,:,2)); 
wy=min(max(wy,1),M);

for x1=1:N,              
 for y1=1:M,
	% Note that wy and wx are opposite to fit with rows and columns!
  warpedim(x1,y1) = im(wy(x1,y1),wx(x1,y1)); 
 end
end
                                                                                                                                                                                                                                                                                                         ./bugger.m                                                                                          000644  017223  001747  00000000036 07365035632 014177  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function bugger;

disp('Yep')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ./ca.m                                                                                              000644  002206  001747  00000000012 07462034005 013463  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         close all
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ./calc_rms_diff.m                                                                                   000664  002226  001747  00000001510 07323360060 015612  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [rms,angl,transl] = calc_rms_diff(tr1,tr2,max_r,centre)
%  [rms,angl,transl] = CALC_RMS_DIFF(tr1,tr2,max_r,centre)
%
%  Calculates the mean squared error between the two candidate
%   transformations (4x4 affine matrices)
%  Returns the rms error (rms), the angle difference (angl) and
%   the translational difference (transl)
%  The domains can be voxel domains (must be the same for both tr1 and tr2)
%   but the range must be measurement space (i.e. in mm)
%  Must specify the maximum radius (max_r), the centre of rotation (world coords)
%   and the measurement sampling (sr - a diagonal matrix)
%

[r,t,s,k]=decomp_aff(tr1*inv(tr2),centre);
angl=real(decomp_rot(r));
isodifft = tr1*inv(tr2) - eye(4);
dn=isodifft(1:3,1:3);
transl=isodifft(1:3,4) + dn*centre;
rms = sqrt( dot(transl,transl) + 1/5*max_r^2 * trace(dn.'*dn) );

                                                                                                                                                                                        ./clamp.m                                                                                           000664  002226  001747  00000000563 07323360060 014142  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function cdata = clamp(data,imin,imax)
%  cdata = clamp(data,imin,imax)
% 
% Returns a clamped copy (cdata) of the input data (data) where
%  each value is clamped to lie between imin and imax

if (imax<imin),
  error('imax must be greater than or equal to imin');
end

mask0 = data<imin;
mask1 = data>imax;
cdata = (1-mask1).*(1-mask0).*data + imin*mask0 + imax*mask1;

                                                                                                                                             ./combine_eros_xfms.m                                                                               000664  002226  001747  00000001765 07323360060 016554  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         % Combines the various transformation files associated with the EROS
%  corrected AIR transforms
%
% It takes the resampling from E0.xfm, the eros transform from E.xfm
%  and the AIR transform from E2.xfm
% The output transformation is written as a shadow transformation in
%  Eros.xfm

cd /usr/people/steve/reg/stroke_stpendle
dirs=str2mat('ab980121sk.Ly1','ag050398.Mf1','ag120897sk.Mm2');
dirs=str2mat(dirs,'ag220498sk.N11','at980108sk.Ll1','bd900113s.Lq1');
dirs=str2mat(dirs,'cd030498.MI1','dg230697sk.I71','dw160697sk.I01');
dirs=str2mat(dirs,'dy190697sk.I31','jm980212_s.LU1','kd040398sk.Me1');
dirs=str2mat(dirs,'mc120697sk.HW1','nb300697sk.Ie1','ps190298sk.M11');
dirs=str2mat(dirs,'rb120697sk.HW1','sf040298sk.MH1','wp260298sk.M81');
for n=1:length(dirs(:,1)),
  eval(['cd ',dirs(n,:)]);
  [t1,dum1,dum2]=read_medx_xfm('E0.xfm');
  [t2,dum1,dum2]=read_medx_xfm('E.xfm');
  [dum1,t3,dum2]=read_medx_xfm('E2.xfm');
  t_total=t3*t2*t1;
  write_medx_xfm('Eros.xfm',[172 220 156],[],t_total,[]);
  cd ..
end
           ./confhyp.m                                                                                         000644  017223  001747  00000001035 07441444537 014375  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function z=confhyp(x,a,b,n)
% CONFHYP Computes the confluent hypergeometric function 
% using a series expansion:
%    f(x;a,b)= sum i (G(a)G(b+i)x^i)/(G(a+i)G(b)i),
% where G is the Gamma function.  Notice that f(x;a,a)=exp(x).
% This function solves Kummer's Equation:
%     xf''(x)+(a-x)f'(x)-bf(x)=0,
% with f(0)=0.
% The parameter n should be a scalar that determines how many terms are
% used in the series expansion.
%
% TB 2002
  z=1;
  a=a+n-1; 
  b=b+n-1;
  while n>0
    z=z.*((a./b).*(x/n))+1;
    a=a-1; b=b-1; n=n-1;
  end




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ./cyclic.m                                                                                          000664  002226  001747  00000000315 07323360060 014307  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [beta] = cyclic(alpha,n)
%  [beta] = cyclic(alpha,n)
%  Cyclic performs a cyclic shift on the data by n values

n=-n;
m=length(alpha(1,:));
if n<0,
	n=n+m;
end
beta=[alpha(:,n+1:m) alpha(:,1:n)];
                                                                                                                                                                                                                                                                                                                   ./cyclic2.m                                                                                         000664  002226  001747  00000000262 07323360060 014372  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [beta] = cyclic2(alpha,p,q)
% function [beta] = cyclic2(alpha,p,q)
% Cyclic performs a cyclic shift of (p,q) on the data

beta=cyclic(alpha,q);
beta=cyclic(beta.',p).';
                                                                                                                                                                                                                                                                                                                                              ./decomp_aff.m                                                                                      000664  002226  001747  00000002310 07323360060 015121  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [rotmat,transl,scales,skew] = decomp_aff(transmat,varargin)
% [rotmat,transl,scales,skew] = DECOMP_AFF(transmat,centre)
%
%  Takes the 4x4 transformation matrix and returns the 3x3 rotation 
%   matrix (rotmat), the 3x1 translation matrix (transl),
%   and the 4x4 matrixes for scaling (scales) and skew (skew)
%   such that:  transmat = [rotmat transl; 0 0 0 1]*skew*scales
%  Centre is an optional argument specifying the centre of the
%   rotation - this affects the translation  (default centre is [0 0 0])
%   with a non-zero centre : 
%      transmat = [rotmat transl + (eye(3)-rotmat)*centre; 0 0 0 1]*skew*scales
%
if (length(varargin)<1)
  centre=[0 0 0].';
else
  centre=varargin{1};
  csz=size(centre);
  if (csz(1)<csz(2)),
    centre=centre.';
  end
end
transl=transmat(1:3,4);
affmat=transmat(1:3,1:3);
x=affmat*[1;0;0];
y=affmat*[0;1;0];
z=affmat*[0;0;1];
sx=norm(x);
sy=sqrt(dot(y,y) - dot(x,y)^2/(sx^2));
a = dot(x,y)/(sx*sy);
x0=x/sx;
y0=y/sy - a*x0;
sz=sqrt(dot(z,z) - dot(x0,z)^2 - dot(y0,z)^2);
b = dot(x0,z)/sz;
c = dot(y0,z)/sz;
scales=diag([sx sy sz 1]);
skew=[1 a b 0; 0 1 c 0; 0 0 1 0; 0 0 0 1];
rotmat=affmat*inv(scales(1:3,1:3))*inv(skew(1:3,1:3));
transl=transl-(eye(3)-rotmat)*centre;
                                                                                                                                                                                                                                                                                                                        ./decomp_rot.m                                                                                      000664  002226  001747  00000001045 07323360060 015175  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [angles] = decomp_rot(rotmat)
%  [angles] = DECOMP_ROT(rotmat)
%
%  Takes the 3x3 rotation matrix and returns the 
%   axis * angle (in radians)
%
rot=rotmat(1:3,1:3);
residual=sum(sum((rot*rot.' - eye(3)).^2));
residual=residual + sum(sum((rot.'*rot - eye(3)).^2));
if (residual>1e-4)
  disp(['WARNING: Failed orthogonality check (residual = ',num2str(residual),')']);
  disp('         Matrix is not strictly a rotation matrix');
  %return;
end
[v,d]=eig(rot);
[y,idx]=sort(real(diag(d)));
angles=abs(angle(d(idx(2),idx(2))))*v(:,idx(3));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ./digamma.m                                                                                         000644  017223  001747  00000000652 07417024067 014325  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function dg=digamma(x,h)
% DIGAMMA - numerical approximation to digamma function 0.
% 
% DG=DIGAMMA(X)
% DG=DIGAMMA(X,H)
% 
% evaluates gamma'(x)/gamma(x) using a numerical approximation 
% to the derivative. Step size in the numerical approximation is H
% (default 0.00001)
%
%        (C) T. Behrens 2002 

if(nargin==1);h=0.00001;end
gamdash=(gamma(x-2*h)-8*gamma(x-h)+8*gamma(x+h)-gamma(x+2*h))/12/h;
dg=gamdash./gamma(x);
                                                                                      ./digammaint.m                                                                                      000644  017223  001747  00000000572 07417051732 015040  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function dg=digammaint(x,h)
% DIGAMMAINT - calculates digamma 0 of an integer
% 
% DG=DIGAMMA(X)
% DG=DIGAMMA(X,H)
% 
% evaluates gamma'(x)/gamma(x) using a numerical approximation 
% to the derivative. Step size in the numerical approximation is H
% (default 0.00001)
%
%        (C) T. Behrens 2002 
x=round(x);
for(i=1:length(x))
  dg(1,i) = digamma(1)+sum(1./(1:x(i)-1));
end                                                                                                                                      ./dirichlet.m                                                                                       000644  017223  001747  00000001444 07417620474 014701  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function pdf=dirichlet(x,u);
% DIRICHLET - calculate the dirichlet pdf
%
%    PDF=DIRICHLET(X,U);
%
%    X is vector of values. 0<X(i)<1 and
%    sum(X,1) = 1
%
%    U is vector of dirichlet parameters
%    0<U(i) and size(U)=[size(X,1) 1];
%
%    dirichlet pdf is ::
%
%    gamma(sum(U,1))/prod(gamma(U),1)*prod(X.^(U-1));
%    for each column of X 
%    (C) T.Behrens 2001
if(~isempty(find(x<=0)) | ~isempty(find(x>=1)))
  error(' X is out of range ')
elseif(~isempty(find(u<=0)))
  error(' U may not be negative or zero ')
elseif(size(u,1)~=size(x,1))
  error(' U must have the same no of rows as X ')
elseif(size(u,2)~=1)
  error(' U may only have one column - define different parameters in different rows of U')
end

u=repmat(u,1,size(x,2));
pdf=gamma(sum(u,1))./prod(gamma(u),1).*prod(x.^(u-1));

                                                                                                                                                                                                                            ./dirichlet2.m                                                                                      000644  017223  001747  00000001465 07417047435 014766  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function pdf=dirichlet2(x,u);
% DIRICHLET2 - calculate the dirichlet pdf
%
%    PDF=DIRICHLET2(X,U);
%
%    X is vector of values. 0<X(i)<1 and
%    sum(X,1) = 1
%
%    U is vector of dirichlet parameters
%    0<U(i) and size(U)=[size(X,1) 1];
%
%    dirichlet pdf is ::
%
%    gamma(sum(U,1))/prod(gamma(U),1)*prod(X.^(U-1));
%    for each column of X 
%    (C) T.Behrens 2001
if(~isempty(find(x<=0)) | ~isempty(find(x>=1)))
  error(' X is out of range ')
elseif(~isempty(find(u<=0)))
  error(' U may not be negative or zero ')
elseif(size(u,1)~=size(x,1))
  error(' U must have the same no of rows as X ')
elseif(size(u,2)~=1)
  error(' U may only have one column - define different parameters in different rows of U')
end

u=repmat(u,1,size(x,2));
pdf=exp(gammaln(sum(u,1)) - sum(gammaln(u),1) + sum(log(x.^(u-1))));

                                                                                                                                                                                                           ./dq.m                                                                                              000644  002206  001747  00000000007 07462034034 013512  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         dbquit
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         ./edim.m                                                                                            000664  002226  001747  00000000572 07323360060 013764  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [d] = edim(a)
% EDIM   Returns the effective dimension of a variable
%    EDIM(A) returns the effective dimension of the variable A
%    That is, if A is 5x4x2 then it returns 3
%    However, if A is 5x1x2 then it returns 2
%    Similary, if A is 2x1 it returns 1 and if A is 1x1 it returns 0
%
%    See also DIM
nulldims = sum(size(a)==1);
d=length(size(a))-nulldims;
                                                                                                                                      ./evsave.m                                                                                          000644  002203  001747  00000000465 07516015506 013643  0                                                                                                    ustar 00heidi                           analysis                        000000  000000                                                                                                                                                                         function evsave(fname,ev);
% EVSAVE saves ev to text.
% 
% EVSAVE(fname,ev);
% 
% FNAME is filename 
%
% EV is ev variable
%
% HJB
if(nargin<2)
  error('Not enough input arguments.')
end

fid=fopen(fname,'w');
if(fid==-1)
  error('Cannot write file')
else
  fprintf(fid,'%12.8f\n',ev);
  fclose(fid);
end

  
                                                                                                                                                                                                           ./fftseries.m                                                                                       000664  002226  001747  00000000275 07323360060 015040  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function fftseries(tseries, i, j, k, varargin)

tseries = detrend(tseries,0);
ft = abs(fft(tseries));
plot(ft(1:ceil(length(tseries)/2)));

title(sprintf('Slice %d; X: %d  Y: %d',k,i,j)); 
                                                                                                                                                                                                                                                                                                                                   ./flipy.m                                                                                           000644  017223  001747  00000000523 07370007005 014035  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function B=flipy(A);
%FLIPY Flip 2D,3D,4D object along y direction
%
%  B=FLIPY(A)
%
%  Tim Behrens

I1=1:size(A,2);
I2=size(A,2):-1:1;
B=zeros(size(A));
if(ndims(A)==2)
  B(:,I1)=A(:,I2);
elseif(ndims(A)==3)
  B(:,I1,:)=A(:,I2,:);
elseif(ndims(A)==4)
  B(:,I1,:,:)=A(:,I2,:,:);
else
  error('Sorry - only works for 2D,3D,4D objects')
end
                                                                                                                                                                             ./gammapdf.m                                                                                        000644  017223  001747  00000000607 07434756326 014513  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function val = gammapdf(a, b, x)

% function val = gammapdf(a, b, x)
%
% Evaluates a gamma of parameters a,b for all values in vector x
% Will have mean m and variance v: m = a/b; v = a/b^2;
% so, a = m^2/v; b = m/v;

% equation is   pow(l,h) * pow(x,h-1) * exp(-l*x - gammln)  but
% replace with more numerically sensible equivalent.

val = exp(a*log(b) + (a-1)*log(x) - b*x - gammaln(a));
                                                                                                                         ./gauss.m                                                                                           000775  002226  001747  00000001665 07323360060 014177  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [g] = gauss(sigma,varargin)
% GAUSS   Creates a gaussian function
%    GAUSS(sigma,n) is an n-by-n gaussian matrix (st.dev = sigma)
%    GAUSS(sigma,m,n) is an m-by-n gaussian matrix (st.dev = sigma)
%    The result is normalised so that the Gaussian sums to 1

if (length(varargin)<1),
  disp('GAUSS must have at least two parameters');
  return;
end
m=varargin{1};
if (length(varargin)>1),
  n=varargin{2};
else
  n=m;
end
mx=(n+1)/2;
my=(m+1)/2;
x=kron(ones(m,1),(1:n)-mx);
y=kron(ones(1,n),(1:m).'-my);
if (length(varargin)<=2),
  g=exp(-(x.^2+y.^2)/(2*sigma^2));
  g=g/sum(sum(g));
  return;
end
if (length(varargin)>3),
  disp('GAUSS cannot create Gaussians in more than 3 dimensions... yet')
  return;
end

p=varargin{3};
mz=(p+1)/2;
bigz=zeros(m,n,p);
bigx=zeros(m,n,p);
bigy=zeros(m,n,p);
for q=1:p,
  bigz(:,:,q)=q-mz;
  bigx(:,:,q)=x;
  bigy(:,:,q)=y;
end
g=exp(-(bigx.^2 + bigy.^2 + bigz.^2)/(2*sigma^2));
g=g/sum(sum(sum(g)));
                                                                           ./generate_latent_ar.m                                                                              000644  002206  001747  00000003545 07421071651 016744  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function [x, arcoeff, sarx, sary, sarz] = generate_latent_ar(Nx, Ny, Nz, T, arcoeff, sarx, sary, sarz, sdevs, mn, activation, regressor)
keyboard;

% x_t=s_t+e_t
% s_t=a*s_t-1+n_t
% N=4; T=180; order=3; arcoeff = [0.3, 0.3, -0.3]; [x, ars] = generate_latent_ar(N,N,1,T, arcoeff, 0.1, 100, 1000); save_avw(x,'/usr/people/woolrich/bpm/data/art_ar3_small','s',[4 4 7 3]); save_avw(arcoeff,'/usr/people/woolrich/bpm/data/arcoeff_ar3','s',[4 4 7 3]); 

if(length(size(arcoeff))<3),
  if(size(arcoeff,1)==1)arcoeff=arcoeff';end;
  order=max(size(arcoeff))
  arcoeff=reshape(shiftdim(reshape(repmat(arcoeff,Nx*Ny*Nz,1),order,Nx,Ny,Nz),1),Nx,Ny,Nz,order);
else
  order=size(arcoeff,length(size(arcoeff)));
end;

if(prod(size(sarx))==1), sarx = ones(Nx,Ny,Nz)*sarx; end;
if(prod(size(sary))==1), sary = ones(Nx,Ny,Nz)*sary; end;
if(prod(size(sarz))==1), sarz = ones(Nx,Ny,Nz)*sarz; end;
if(prod(size(sdevs))==1), sdevs = ones(Nx,Ny,Nz)*sdevs; end;
if(prod(size(activation))==1), activation = ones(Nx,Ny,Nz)*activation; end;
if(prod(size(regressor))==1), regressor = zeros(T,1); end;

x = randn(Nx,Ny,Nz,T);
%x = zeros(Nx,Ny,Nz,T);

for t = 1:T,
for i = 1:Nx
for j = 1:Ny,
for k = 1:Nz,

  s(i,j,k,t) = activation(i,j,k)*regressor(t) + randn*sdevs(i,j,k);

  if t > 1,
   if(i>1) s(i,j,k,t) = sarx(i-1,j,k)*s(i-1,j,k,t-1) + s(i,j,k,t); end;
   if(i<Nx) s(i,j,k,t) = sarx(i,j,k)*s(i+1,j,k,t-1) + s(i,j,k,t); end;
   if(j>1) s(i,j,k,t) = sary(i,j-1,k)*s(i,j-1,k,t-1) + s(i,j,k,t); end;
   if(j<Ny) s(i,j,k,t) = sary(i,j,k)*s(i,j+1,k,t-1) + s(i,j,k,t); end;
   if(k>1) s(i,j,k,t) = sarz(i,j,k-1)*s(i,j,k-1,t-1) + s(i,j,k,t); end;
   if(k<Nz) s(i,j,k,t) = sarz(i,j,k)*s(i,j,k+1,t-1) + s(i,j,k,t); end;

   for a = 1 : order,
	if t > a,
           s(i,j,k,t) = arcoeff(i,j,k,a)*s(i,j,k,t-a) + s(i,j,k,t);
           x(i,j,k,t) = s(i,j,k,t) + x(i,j,k,t);
	end;
   end;
  end;
end;
end;
end;
end;

x = x+mn;	

                                                                                                                                                           ./generate_spm_xfms.m                                                                               000664  002226  001747  00000001470 07323360060 016552  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         % Generates a t2_spm_rot0.xfm file for each directory, by
%  converting the t2_spm_rot0.m file into medx format

init_size=[256 256 30];
final_size=[172 220 156];
cd /usr/people/steve/reg/stroke_stpendle
dirs=str2mat('ab980121sk.Ly1','ag050398.Mf1','ag120897sk.Mm2');
dirs=str2mat(dirs,'ag220498sk.N11','at980108sk.Ll1','bd900113s.Lq1');
dirs=str2mat(dirs,'cd030498.MI1','dg230697sk.I71','dw160697sk.I01');
dirs=str2mat(dirs,'dy190697sk.I31','jm980212_s.LU1','kd040398sk.Me1');
dirs=str2mat(dirs,'mc120697sk.HW1','nb300697sk.Ie1','ps190298sk.M11');
dirs=str2mat(dirs,'rb120697sk.HW1','sf040298sk.MH1','wp260298sk.M81');
for n=1:length(dirs(:,1)),
  eval(['cd ',dirs(n,:)]);
  sxfm=get_trans('t2_spm_rot0.m');
  mxfm=spm2medx(sxfm,init_size,final_size);
  write_medx_xfm('t2_spm_rot0.xfm',final_size,[],mxfm,[]);
  cd ..
end
                                                                                                                                                                                                        ./genrot.m                                                                                          000644  017223  001747  00000000615 07377710030 014220  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function R=genrot(V,th);
%GENROT - generic rotation matrix - Matrix for 3D rotation about V
%through angle th.
%       R=GENROT(V,TH)
%
% Tim Behrens

if nargin~=2
  error('fool')
end
V = V/sqrt(sum(V.^2));
C = cos(th); S = sin(th);
x=V(1);y=V(2);z=V(3);

R=[x^2+C*(1-x^2) x*y*(1-C)-z*S z*x*(1-C)+y*S;
   x*y*(1-C)+z*S y^2+C*(1-y^2) y*z*(1-C)-x*S;
   z*x*(1-C)-y*S y*z*(1-C)+x*S z^2+C*(1-z^2)];


                                                                                                                   ./get_trans.m                                                                                       000664  002226  001747  00000000461 07323360060 015031  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [transmat] = get_trans(fname)
% [transmat] = get_trans(fname)
% Opens the file fname and returns the 4x4 transformation matrix
fid = fopen(fname);
rmat = fscanf(fid,'%f');
fclose(fid);
if (prod(size(rmat))==16),
  transmat=reshape(rmat,4,4).';
else
  transmat=[reshape(rmat,4,3).'; 0 0 0 1];
end
                                                                                                                                                                                                               ./getlboxvals.m                                                                                     000664  002226  001747  00000001476 07323360060 015404  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function getlboxvals(vol,num)
% getlboxvals(vol,num) allows num individual coordinates and intensities to be
%  found from a lightbox plot

sx = length(vol(:,1,1));
sy = length(vol(1,:,1));
sz = length(vol(1,1,:));
sz1=round(sqrt(sz));
sz2 = ceil(sz/sz1);

if (num>=0),
  count=1;
else
  count=num-1;
end

button = 1;
while ((button<=1) & (count<=num)),
  [xp,yp,button] = ginput(1);
  xc = round(xp);
  yc = round(yp);
  axh = gca;

  if ( (xc>0) & (yc>0) & (xc<=sx) & (yc<=sy) ),
    for y=1:sz1,
      for x=1:sz2,
	imno = x+(y-1)*sz2;
	subplot(sz1,sz2,imno);
	if (imno<=sz),
	  if (axh==gca),
	    disp(['X Y Z ; I = ',num2str(round(yp)),' ', ...
		  num2str(round(xp)),' ', ...
		  num2str(imno),' ; ',num2str(vol(round(yp),round(xp),imno))]);
	  end
	end
      end
    end
  end
  if (num>=0),
    count=count+1;
  end
end

                                                                                                                                                                                                  ./glm_fit.m                                                                                         000664  002206  001747  00000002224 07516035120 014531  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function glm_fit(tseries, i, j, k, copenum, t_or_f)

% glm_fit(tseries, i, j, k)
%
% used by ViewFMRI in view_feat_results

global xs dm copes zs acs c fc;

cla;
hold on;

filter = establish_filter(squeeze(acs(i,j,k,:)),length(tseries));

if(std(filter)>0),
ts = prewhiten(detrend(tseries,0),filter);
plot(ts);
pwdm = prewhiten(dm,filter);

if(t_or_f=='f')
  j=1;
  for i=1:size(fc,2),
    if(fc(copenum,i))
      c1(j,:) = c(i,:);
      j=j+1;
    end;
  end;
else
  c1 = c(copenum,:);
end;

% get effective regressor(s) for contrast(s):
ev = (c1*pinv(pwdm))';

b = pinv(ev)*ts;
plot(ev*b,'r');
end;

title(sprintf('Slice %d; X: %d  Y: %d Zstat:4.2f',i,j,k,zs(i,j,k))); 

hold off;

%%%%%%%%%%

function y = prewhiten(x,filter)

n = size(x,1);
y = x;

for i = 1:size(x,2),
	zeropad = length(filter);
	fx = fft(squeeze(x(:,i)),zeropad);
	y2 = real(ifft(fx.*filter));
	y(:,i)=y2(1:n);
end;

%%%%%%%%%%%

function fy = establish_filter(ac,n)

zeropad = 2^nextpow2(n);
ac1 = zeros(zeropad,1);
ac1(1:length(ac)) = ac;
ac1(zeropad-length(ac)+2:zeropad) = fliplr(ac(2:length(ac)));

fac = fft(ac1,zeropad);

warning off;
fac = 1./fac;
warning on;

fac(1) = 0;
fy = fac/std(fac);                                                                                                                                                                                                                                                                                                                                                                            ./glm_fit2.m                                                                                        000644  002206  001747  00000001764 07360325771 014634  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function glm_fit2(tseries, i, j, k, copenum)

% glm_fit(tseries, i, j, k)
%
% used by ViewFMRI in view_feat_results

global xs dm pes zs acs c;

cla;
hold on;

filter = establish_filter(squeeze(acs(i,j,k,:)),length(tseries));

if(std(filter)>0),
ts = prewhiten(detrend(tseries,0),filter);
%plot(ts);
pwdm = prewhiten(dm,filter);

ev = (c(1,:)*pinv(pwdm))';
for cn=2:4,
 ev = ev + (c(cn,:)*pinv(pwdm))';
end;

plot(real(fft(ts-ev)),'r')
end;

title(sprintf('Slice %d; X: %d  Y: %d Zstat:4.2f',i,j,k,zs(i,j,k))); 

hold off;

%%%%%%%%%%

function y = prewhiten(x,filter)

n = size(x,1);
y = x;

for i = 1:size(x,2),
	zeropad = length(filter);
	fx = fft(squeeze(x(:,i)),zeropad);
	y2 = real(ifft(fx.*filter));
	y(:,i)=y2(1:n);
end;

%%%%%%%%%%%

function fy = establish_filter(ac,n)

zeropad = 2^nextpow2(n);
ac1 = zeros(zeropad,1);
ac1(1:length(ac)) = ac;
ac1(zeropad-length(ac)+2:zeropad) = fliplr(ac(2:length(ac)));

fac = fft(ac1,zeropad);

warning off;
fac = 1./fac;
warning on;

fac(1) = 0;
fy = fac/std(fac);            ./gmmseries.m                                                                                       000644  002206  001747  00000000503 07441454362 015111  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function gmmseries(tseries, i, j, k, varargin)

if(nargin>4)
n=varargin;
else
n=1;
end;

if std(tseries)<=0,return;end;

if(n>1)
fit_gmm2(tseries,n,1);
else,
mix = gmm(1, 1, 'spherical');
mix.centres = mean(tseries);
mix.covars = var(tseries);
plot_mm(mix,tseries,1);
end;

title(sprintf('Slice %d; X: %d  Y: %d',k,i,j)); 
                                                                                                                                                                                             ./hello.m                                                                                           000644  017223  001747  00000000056 07377477107 014043  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function hello;
disp('I am a computer. Fool.')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ./hist2d.m                                                                                          000664  017223  001747  00000003276 07464462276 014144  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function [bin,vecx,vecy] = hist2d(mat,x,y);

% bin = hist2d(mat);
% 2-D histogram with 10*10 bins
% mat should be N*2, where first col is x and second is y
%
% bin = hist2d(mat,n);
% 2-D histogram with n*n bins
%
% bin = hist2d(mat,nx,ny);
% 2-D histogram with nx*ny bins
%
% bin = hist2d(mat,x);
% 2-D histogram with bins specified by x in both directions 
% (x must be regularly spaced)
%
% bin = hist2d(mat,x,y);
% 2-D histogram with bins specified by x and y 
% (x and y must be regularly spaced)
%
% Copyright(C) 2001
%  M.Jenkinson, T.Behrens
%
% The authors acknowledge that there was a very small contribution
% from one other member of the FMRIB image analysis group, but shy away
% from mentioning his name to avoid causing unnecessary angst to
% users who have previously had trouble with the countless bugs in
% ALL of his other code.

if(size(mat,2) ~=2), mat=mat'; end;

if(nargin == 1),
  x=10;
end;

if(length(x) == 1)
  nx=x;
  
  if(nargin <= 2),
    ny=x;
  else,
    ny=y;
  end;

  maxx=max(mat(:,1));minx=min(mat(:,1));rx=maxx-minx;
  maxy=max(mat(:,2));miny=min(mat(:,2));ry=maxy-miny;
else,
  nx = length(x);	
  maxx=max(x);minx=min(x);rx=maxx-minx;
  if(nargin == 2),
    ny=length(x);
    maxy=max(x);miny=min(x);ry=maxy-miny;
  else,
    ny=length(y);	
    maxy=max(y);miny=min(y);ry=maxy-miny;
  end;
end;

bin = zeros(nx,ny);

for i=1:length(mat);
  xbin=ceil(nx*(mat(i,1)-minx)/rx);
  ybin=ceil(ny*(mat(i,2)-miny)/ry);
  if(~(xbin<1 | ybin<1 | xbin>nx | ybin>ny))
    bin(xbin,ybin)=bin(xbin,ybin)+1;
  end;
end
vecy=[miny:ry/(ny-1):maxy];
vecx=[minx:rx/(nx-1):maxx];
%keyboard
if(nargout==0)
  imagesc(vecy,vecx,bin(1:size(bin,1),1:size(bin,2)));
  colorbar;
  axis xy;axis image;
end










                                                                                                                                                                                                                                                                                                                                  ./hist2dsph.m                                                                                       000644  017223  001747  00000001174 07464273364 014646  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function [x,y,z]=hist2dsph(thph,n);
%HIST2DSPH -- 2D histogram of theta, phi round a sphere
%
%        [X,Y,Z]=HIST2DSPH(THPH,N);
%
%        THPH is Mx2 matrix of Theta and Phi
%
%        N is number of bins on each axis
%       
% Copyright (C) T.Behrens 2001


bin=hist2d(thph,n);
th=min(thph(:,1)):(max(thph(:,1))-min(thph(:,1)))/(n-1):max(thph(:,1));
ph=min(thph(:,2)):(max(thph(:,2))-min(thph(:,2)))/(n-1):max(thph(:,2));
x=sin(th')*cos(ph); y=sin(th')*sin(ph); z=cos(th')*ones(1,n);
x=x.*bin; y=y.*bin; z=z.*bin;
if nargout==0
  h=surf(x,y,z); set(h,'linestyle','none'); colormap([1 0 1]);
  light; lighting phong; axis equal;
end                                                                                                                                                                                                                                                                                                                                                                                                    ./histseries.m                                                                                      000644  002226  001747  00000000210 07354343151 015221  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function histseries(tseries, i, j, k, varargin)

hist(tseries,sqrt(length(tseries)));

title(sprintf('Slice %d; X: %d  Y: %d',k,i,j)); 
                                                                                                                                                                                                                                                                                                                                                                                        ./ho.m                                                                                              000644  002206  001747  00000000010 07461744555 013525  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         hold on;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ./hoff.m                                                                                            000644  017223  001747  00000000012 07462251776 013647  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         hold off;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ./how.m                                                                                             000644  017223  001747  00000000044 07352653215 013517  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function how;
disp('Ask Jenkybaby');                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ./imagesc0.m                                                                                        000775  002226  001747  00000001125 07323360060 014534  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function imagesc0(varargin)
% IMAGESC0(...,im)
%
%   Displays the image (im) using the IMAGE command
%     but with a pre-scaling to scale the intensity
%     to fit the maximum range in the colormap with
%     the value zero ALWAYS being the middle of the
%     colormap
%   All initial arguements are passed directly to IMAGE

n = length(varargin);
if (n==0),
   error('Imagesc must have at least one arguement');
end
im = varargin{n};
mx=max(max(im));
mn=min(min(im));
me=mean(mean(im));
sc=max(abs(mx-me),abs(mn-me));
image(varargin{1:(n-1)},((im-me)/sc + 1)*0.5*max(size(colormap)))
axis off
                                                                                                                                                                                                                                                                                                                                                                                                                                           ./imagesc1.m                                                                                        000775  002226  001747  00000000713 07323360060 014537  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function imagesc(varargin)
% IMAGESC(...,im)
%
%   Displays the image (im) using the IMAGE command
%     but with a pre-scaling to scale the intensity
%     to fit the maximum range in the colormap
%   All initial arguements are passed directly to IMAGE

n = length(varargin);
if (n==0),
   error('Imagesc must have at least one arguement');
end
im = varargin{n};
mx=max(max(im));
mn=min(min(im));
image(varargin{1:(n-1)},(im-mn)/(mx-mn)*max(size(colormap)))
                                                     ./imgoverlay.m                                                                                      000644  017223  001747  00000002142 07466461312 015102  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function [cmap,img]=imgoverlay(imga,imgb,cmap1,cmap2);
% IMGOVERLAY - overlay a thresholded image (imgb) on a background (imga)
%
% [cmap,img] = imgoverlay(imga,imgb);
%
% [cmap,img] = imgoverlay(imga,imgb,cmap1);
%
% e.g imgoverlay(imga,imgb,'bone') - background will be bone 
%
% [cmap,img] = imgoverlay(imga,imgb,cmap1,cmap2);
%
% e.g imgoverlay(imga,imgb,'bone','hot') - background bone,
% overlay hot.
% 
% To display - use image or lbox. 
% e.g 
% 
% image(img(:,:,slice_no));colormap(cmap);
%
% lbox(img);colormap(cmap);
%
% orthoview doesn't seem to work - 
% I think MJ does some kind of robust range finding
% 
% TB (2002)
%

if(nargin<2)
  error('Not enough Input arguments');
elseif(~(prod(size(imga)==size(imgb))==1))
  error('Images must be same size -- dickhead');
end

if(nargin<4)
  cmap2='hot';
end
if(nargin<3)
  cmap1='bone';
end
figure;
colormap(cmap1);
c1=colormap;
c1=c1(1:2:end,:);
colormap(cmap2);
c2=colormap;
c2=c2(1:2:end,:);
close;
cmap=[c1;c2];

imgasc=(imga-min3(imga))/range(squash(imga))*32;
imgasc(find(imgb>0))=32;
imgbsc=(imgb-min3(imgb))/range(squash(imgb))*32;
img=imgasc+imgbsc;




                                                                                                                                                                                                                                                                                                                                                                                                                              ./insert.m                                                                                          000644  017223  001747  00000001346 07522257566 014244  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function [MatA]=insert(MatA,MatB,dim,num)
% INSERT - insert rows or columns into a sorted matrix
%
% MATC = INSERT(MATA,MATB,DIM,NUM);
%
% MATA - matrix into which rows or cols are to be inserted
%
% MATB - rows or cols to be inserted
%
% DIM - 1 for rows, 2 for cols
%
% NUM row or column number according to which MATA is sorted
%
% MATC - Output;
%
% TB

if(dim==2)
  MatA=MatA';
  MatB=MatB';
elseif(dim~=1)
  error('Invalid value for dim - must be 1 or 2');
end

for i=1:size(MatB,1)
  I=find(MatA(:,num)<MatB(i,num));
  if(isempty(I))
    MatA=[MatB(i,:);MatA];
  elseif(length(I)==size(MatA,1));
    MatA=[MatA;MatB(i,:)];
  else
  MatA=[MatA(I,:);MatB(i,:);MatA(I(end)+1:end,:)];
  end
end

if(dim==2)
  MatA=MatA';
  MatB=MatB';
end
                                                                                                                                                                                                                                                                                          ./interactivewarp.m                                                                                 000664  002226  001747  00000003171 07323360060 016253  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [warpedim] = interactivewarp(im,smooth)
% [warpedim] = interactivewarp(im,smooth)
%
% Allows the user to specify points, with the mouse, to push and pull
%  the image around
% The optional smoothing factor specifies how smooth the warps will
%  be with 1.0 (the default) being reasonably smooth
% Stops when the return key is pressed

if (nargin<2),
  smooth=1.0;
end

[N,M] = size(im);
warp = zeros(N,M,2);
warp2 = warp;
warpedim = im;
x = kron(1:M,ones(N,1));
y = kron((1:N).',ones(1,M));

imagesc(im);

while 2>1 , 
  undo = 0;

  c=1;
  while c<=2,
   [gx,gy,button]=ginput(1);
   if (length(gx)~=1),
    c=-1;
    break;
   end
 
   if (button==3),
    if (undo==0),
      warp2 = warp;
      undo = 1;
    else
      warp2 = zeros(size(warp));
    end
    warpedim = applywarp(im,warp2,x,y);
    imagesc(warpedim);
    c=0;
   end
 
   if (button==2),
    smooth=input('Enter new smoothing factor: ');
    c=0;
   end
 
   if c>=1,
     mx(c) = gx;  my(c) = gy;
   end
   c=c+1;
  end

  warp = warp2;

  if c==-1,
   break;
  end

  sig = smooth*norm([mx(1)-mx(2) , my(1)-my(2)]);
  sig=max(sig,0.0001);
  damp = exp(-(x-mx(2)).^2/(sig.^2)).*exp(-(y-my(2)).^2/(sig.^2));
  warp2(:,:,1) = warp(:,:,1) + (mx(1)-mx(2))*damp;
  warp2(:,:,2) = warp(:,:,2) + (my(1)-my(2))*damp;

  warpedim = applywarp(im,warp2,x,y);
  imagesc(warpedim);

end

% final smoothing
smoothfilter = ones(3,3);
smoothfilter(2,2) = 2;
smoothfilter = smoothfilter/sum(sum(smoothfilter));
warp(:,:,1) = conv2(warp2(:,:,1),smoothfilter,'same');
warp(:,:,2) = conv2(warp2(:,:,2),smoothfilter,'same');

% Return warpedim
warpedim = applywarp(im,warp,x,y);
imagesc(warpedim);
                                                                                                                                                                                                                                                                                                                                                                                                       ./isany.m                                                                                           000644  017223  001747  00000000246 07526671167 014062  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function A = isany(B,vec);
%% ISANY -- returns B==(any of the elements of vec);
%%
%% A = ISANY(B,VEC);

A=zeros(size(B));
for i=1:length(vec)
  A=A|(B==vec(i));
end
                                                                                                                                                                                                                                                                                                                                                          ./isin.m                                                                                            000644  017223  001747  00000000534 07345651361 013672  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function TF = isin(a,roi);
% ISIN return 1 if a is in roi else 0 
%
% e.g  TF = ISIN([1 2 1],[0 0 0;10 10 10]); returns 1
% but  TF = ISIN([1 20 1],[0 0 0;10 10 10]); returns 0
% works with vectors of any length.
%
% copyright (C) T.Behrens 2001

TF=(repmat(roi(1,:),size(a,1),1)<a&a<repmat(roi(2,:),size(a,1),1));
TF=floor(sum(TF,2)/size(TF,2));                                                                                                                                                                      ./lbox.m                                                                                            000664  002226  001747  00000000173 07330034052 014004  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function lbox(vol,opts)
% Calls lightbox

if (nargin>=2),
  lightbox(vol,opts);
elseif (nargin==1),
  lightbox(vol);
end


                                                                                                                                                                                                                                                                                                                                                                                                     ./lboxc.m                                                                                           000664  002226  001747  00000000357 07330034011 014146  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function lboxc(vol,opts)
% Displays a lightbox view in coronal (given an axial input)
%
% See also: lightbox

if (nargin>=2),
  lbox(flipdim(permute(vol,[3 1 2]),1),opts);
elseif (nargin==1),
  lbox(flipdim(permute(vol,[3 1 2]),1));
end


                                                                                                                                                                                                                                                                                 ./lboxs.m                                                                                           000664  002226  001747  00000000360 07330034027 014167  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function lboxs(vol,opts)
% Displays a lightbox view in sagittal (given an axial input)
%
% See also: lightbox

if (nargin>=2),
  lbox(flipdim(permute(vol,[3 2 1]),1),opts);
elseif (nargin==1),
  lbox(flipdim(permute(vol,[3 2 1]),1));
end


                                                                                                                                                                                                                                                                                ./lightbox.m                                                                                        000664  002210  001747  00000003107 07354661565 014052  0                                                                                                    ustar 00mark                            analysis                        000000  000000                                                                                                                                                                         function lightbox(vol,opts)
% lightbox(vol) displays all slices of the volume vol as a lightbox
%  presentation
% lightbox(vol,num) displays lightbox and allows num points to be 
%  interrogated with the mouse (num=-1 gives unlimited number of
%  points; num=0 is default with no points )
% The interrogation can be terminated by a non-left button press or
%  a keyboard press *within the figure window*
% lightbox(vol,[num imin imax]) specifies the absolute intensity range
%  to use imin <= intensity <= imax.   This is useful for bypassing
%  the default robust range when it is inappropriate.
% lightbox(vol,[num 0]) calculates the absolute intensity range
%     short for lightbox(vol,[num min3(vol) max3(vol)]


if (nargin==1),
  num=0;
else
  num=opts(1);
end 

if ((nargin>1) & (length(opts)>=3)),
  irange = [opts(2) opts(3)];
elseif ((nargin>1) & (length(opts)==2)),
  % Use absolute range
  irange = [min(min3(vol)) max(max3(vol))];
else
  % Use the robust min and max
  irange = percentile(vol,[2 98]);
end

if ( abs(irange(2) - irange(1))==0),
  % Use absolute range if robust is no good
  irange = [min(min3(vol)) max(max3(vol))];
end
if ( abs(irange(2) - irange(1))==0),
  irange = [irange(1) irange(1)+1];
end


ioff = irange(1);
iscale = (length(colormap)-1)/(irange(2) - irange(1));

clf

sz = length(vol(1,1,:));
sz1=round(sqrt(sz));
sz2 = ceil(sz/sz1);


for y=1:sz1,
  for x=1:sz2,
    imno = x+(y-1)*sz2;
    if (imno<=sz),
      subplot(sz1,sz2,imno);
      image(max(min((vol(:,:,imno)-ioff)*iscale + 1,length(colormap)),1));
      axis off
    end
  end
end

getlboxvals(vol,num);

                                                                                                                                                                                                                                                                                                                                                                                                                                                         ./line3.m                                                                                           000644  017223  001747  00000000476 07521721366 013746  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function h=line3(R);
% LINE3 - convenience to save typing LINE(R(:,1),R(:,2),R(:,3));
% 
% LINE3(R);
%
% H = LINE3(R);
%
% performs line above;
%
% TB

if(nargin~=1)
  error('Wrong Number of input arguments - Read the help!!')
end
if(nargout==0)
  line(R(:,1),R(:,2),R(:,3));
else 
  h=line(R(:,1),R(:,2),R(:,3));
end
                                                                                                                                                                                                  ./lms.m                                                                                             000775  002226  001747  00000000267 07323360060 013645  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function ms = lmsfn(p)
% ms = lmsfn(p)
% use global variables xlms and ylms
global xlms ylms
% return the median of (ylms - p1*xlms - p2)^2
ms = median((ylms - p(1)*xlms - p(2)).^2);
                                                                                                                                                                                                                                                                                                                                         ./make_rot.m                                                                                        000664  002226  001747  00000001461 07323360060 014645  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [rotmat] = make_rot(angl,centre)
% [rotmat] = make_rot(angl,centre)
% Creates a 4x4 rotation matrix suitable for combination with the
%  affine transformations returned by MEDx
% angl should be a 3-vector representing the axis of rotation and
%  with length equal to the angle of rotation (in radians)
% centre should be the 3D coordinate of the centre for the rotation

if (norm(angl)==0),
  rotmat=eye(4);
  return;
end
ax=angl/norm(angl);
sz=size(ax);
if (sz(1)==3),
  ax=ax.';
end
x1=ax;
x2=[-ax(2),ax(1),0];
if (norm(x2)==0),
  x2=[1 0 0];
end
x3=cross(x1,x2);
x2=x2/norm(x2);
x3=x3/norm(x3);
tx=[x2; x3; x1];
th=norm(angl);
r=[cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
rot=tx.'*r*tx;
sz=size(centre);
if (sz(2)==3),
  centre=centre.';
end
trans=(eye(3)-rot)*centre;
rotmat=[rot, trans; 0 0 0 1]; 
                                                                                                                                                                                                               ./makefcontrast.m                                                                                   000644  002206  001747  00000000602 07526723446 015765  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function fc = makefcontrast(fccoded, tc, fcopenum)
   
% fc = makefcontrast(fccoded, tc, fcopenum)
%
% makes an f contrast matrix for the f contrast no. fcopenum
% from a t contrast matrix and a fccoded 
% matrix, which informs which t contrasts make up each f contrast

   j=1;
   for i=1:size(fccoded,2),
      if(fccoded(fcopenum,i))
	 fc(j,:) = tc(i,:);
	 j=j+1;
      end;
   end;
                                                                                                                              ./makehist.m                                                                                        000664  002226  001747  00000001157 07323360060 014653  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [hist2] = makehist(im1,im2)
%  hist2 = makehist(im1,im2)
%
%  Makes a 2D joint histogram from images im1 and im2.
%  Uses as many entries as intensity values in each, so for
%   compact histograms preprocess the images so that only bin numbers are
%   used
%  Assumes both images have the same dimensions

[M,N]=size(im1);

min1 = floor(min(min(im1))) - 1;
min2 = floor(min(min(im2))) - 1;
hist2 = zeros(ceil(max(max(im1)))-min1, ceil(max(max(im2)))-min2);
for x1=1:M,
 for y1=1:N,
  hist2(round(im1(x1,y1))-min1,round(im2(x1,y1))-min2) = ...
    hist2(round(im1(x1,y1))-min1,round(im2(x1,y1))-min2) + 1;
 end
end
                                                                                                                                                                                                                                                                                                                                                                                                                 ./mandelbrot.m                                                                                      000644  017223  001747  00000001723 07434503474 015060  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function varargout=mandelbrot(varargin)
% MANDELBROT - Calculate mandelbrot set.
% 
%  MANDELBROT(A) - display the set. A is a measure of resolution
%  (default 500).
% 
%  [C] = MANDELBROT(A) - C is output.
%
% TB 2002


if(nargin>2)
  error('Too many input arguments to MANDELBROT')
elseif(nargout>1)
  error('Too many output arguments to MANDELBROT')
end

if nargin==0
  a=500;
elseif nargin ==1
  a=varargin{1};
  plane='mu';
else
  a=varargin{1};
  plane=varargin{2};
end

A=zeros(2*a+1,3*a+1);
C=zeros(2*a+1,3*a+1);
creal=repmat((-2*a:a)/a,2*a+1,1);
cimag=repmat(sqrt(-1)/a*[-a:a]',1,3*a+1);
c=creal+cimag;
for k=1:32
  if(strcmp(plane,'mu'))
    A=A.^2+c;
  elseif(strcmp(plane,'invmu'))
    warning off;
    A=A.^2+1./c;
  elseif(strcmp(plane,'invmu2'));
    warning off
    A=A.^2-1./(c+0.25);
  else
    error('No such plane');
  end
  
    B=abs(A)>4&C==0;
  C=C+k*B;
end

if(nargout==0)
  image([-2*a:a]/a,[-a:a]/a,C);axis image;axis xy;
else
  varargout{1}=C;
end


                                             ./max3.m                                                                                            000664  002226  001747  00000000111 07323360060 013703  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [m] = max3(x)
% performs max(max(max(x)))

m = max(max(max(x)));                                                                                                                                                                                                                                                                                                                                                                                                                                                       ./mean3.m                                                                                           000664  002226  001747  00000000120 07323360060 014036  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [m] = mean3(x)
% performs mean(mean(mean(x)))

m = mean(mean(mean(x)));                                                                                                                                                                                                                                                                                                                                                                                                                                                ./meannz.m                                                                                          000664  002210  001747  00000000315 07326323137 013505  0                                                                                                    ustar 00mark                            analysis                        000000  000000                                                                                                                                                                         function [m] = meannz(x)
% finds the mean of all the non-zero entries
% will return the overall mean regardless of object dimension (uses reshape)

nx=reshape(x,[prod(size(x)) 1]);
m = sum(nx)/sum(nx~=0);
                                                                                                                                                                                                                                                                                                                   ./median3.m                                                                                         000664  002226  001747  00000000162 07323360060 014361  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [m] = median3(x)
% performs median(reshape(x,prod(size(x)),1));

m = median(reshape(x,prod(size(x)),1));
                                                                                                                                                                                                                                                                                                                                                                                                              ./mediannz.m                                                                                        000664  002210  001747  00000000605 07326323542 014024  0                                                                                                    ustar 00mark                            analysis                        000000  000000                                                                                                                                                                         function [m] = mediannz(x)
% finds the median of all the non-zero entries
% will return the overall median regardless of object dimension (uses reshape)

nx=reshape(x,[prod(size(x)) 1]);
p = 100*(1 - 0.5*sum(nx~=0)/length(nx));
% Now set all zero elements to be less than the minimum value, and hence
%  at the start of the percentiles
nx = nx + (min(nx)-1)*(nx==0);
m = percentile(nx,p);
                                                                                                                           ./min3.m                                                                                            000664  002226  001747  00000000111 07323360060 013701  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [m] = min3(x)
% performs min(min(min(x)))

m = min(min(min(x)));                                                                                                                                                                                                                                                                                                                                                                                                                                                       ./mri2medx.m                                                                                        000664  002226  001747  00000002237 07323360061 014576  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [transmat] = mri2medx(tr_mat,sr,init_size,final_size,finalcentre)
% [transmat] = mri2medx(tr_mat,sr,init_size,final_size,finalcentre)
%
% Converts a transformation matrix between mri format
%   and medx format  (pre and post flipped in y and adjusted for centre)
% It must be given: 
%       resampling matrix (sr)  - a diagonal matrix
%       the sizes of the initial and final
%          volumes init_size, final_size both as [xmax,ymax,zmax].
%             Convention is that xmax = no x voxels , ymax = etc.
%       centre of the final volume (finalcentre) - a 3x1 matrix

csze=size(finalcentre);
if (csze(1)<csze(2)),
  finalcentre=finalcentre(1,1:3).';
end
% posttr transforms from world -> minc voxels of the av305
posttr=[diag([1 1 1]), finalcentre; 0 0 0 1];
% pretr transforms from minc voxels -> world
pretr=tr_mat*sr;

flipy=diag([-1 -1 1 1]);
% flip1 transforms from medx voxels -> minc voxels (initial volume)
flip1=flipy;
flip1(1,4)=init_size(1)-1;
flip1(2,4)=init_size(2)-1;
% flip2 transforms from minc voxels -> medx voxels (final volume)
flip2=diag([-1 -1 1 1]);
flip2(1,4)=final_size(1)-1;
flip2(2,4)=final_size(2)-1;

transmat=flip2*posttr*pretr*flip1;
                                                                                                                                                                                                                                                                                                                                                                 ./mvtpdf.m                                                                                          000644  002206  001747  00000000756 07522322015 014415  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function p = mvtpdf(t, m, C, v)

% p = mvtpdf(t, m, C, v)
%
% Type 1 multivariate t pdf
% t is vector for which the probability is to be calculated
% m is the mean vector
% C is the covariance matrix
% v is the dof

if(size(t,1) ~= length(C))
  t = t';
end;
if(size(m,1) == 1)
  m = m';
end;

k = length(m);

ts = t;

for i = 1:size(ts,2),
term = exp(gammaln((v + k) / 2) - gammaln(v/2));
p(i) = term ./ (sqrt(det(C))*(v*pi).^(k/2) .* (1 + (ts(i)-m)'*inv(C)*(ts(i)-m)/ v) .^ ((v + k)/2));
end;
                  ./mydim.m                                                                                           000664  002226  001747  00000000375 07323360061 014167  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [d] = dim(a)
% DIM   Returns the dimension of a variable
%    DIM(A) returns the dimension of the variable A
%    Note that this is the dimension of the MATLAB variable, not
%    its effective dimension.
%
%    See also EDIM
d=length(size(a));
                                                                                                                                                                                                                                                                   ./oldmri2medx.m                                                                                     000664  002226  001747  00000002406 07323360061 015273  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [transmat] = mri2medx(tr_mat,sr,init_size,final_size,finalcentre)
% [transmat] = mri2medx(tr_mat,sr,init_size,final_size,finalcentre)
%
% Converts a transformation matrix between mri format
%   and medx format  (pre and post flipped in y and adjusted for centre)
% It must be given: 
%       resampling matrix (sr)  - a diagonal matrix
%       the sizes of the initial and final
%          volumes init_size, final_size both as [xmax,ymax,zmax].
%             Convention is that xmax = no x voxels , ymax = etc.
%       centre of the final volume (finalcentre) - a 3x1 matrix

csze=size(finalcentre);
if (csze(1)<csze(2)),
  finalcentre=finalcentre(1,1:3).';
end
% posttr transforms from world -> minc voxels of the av305
%xoffset=[final_size(1)-1-finalcentre(1); finalcentre(2:3)];
%posttr=[diag([-1 1 1]), xoffset; 0 0 0 1];
posttr=[diag([1 1 1]), finalcentre; 0 0 0 1];
% pretr transforms from minc voxels -> world
pretr=tr_mat*sr;

flipy=diag([1 -1 1 1]);
% flip1 transforms from medx voxels -> minc voxels (initial volume)
flip1=flipy;
flip1(2,4)=init_size(2)-1;
% flip2 transforms from minc voxels -> medx voxels (final volume)
%flip2=diag([-1 -1 1 1]);
%flip2(1,4)=final_size(1)-1;
flip2=diag([1 -1 1 1]);
flip2(2,4)=final_size(2)-1;

transmat=flip2*posttr*pretr*flip1;
                                                                                                                                                                                                                                                          ./orthoview.m                                                                                       000664  002210  001747  00000003557 07353605471 014262  0                                                                                                    ustar 00mark                            analysis                        000000  000000                                                                                                                                                                         function orthoview(vol,centre,voxdim,opts)
% orthoview(vol) displays orthogonal views (sagittal, coronal and axial) of
%  a volume, using the centre of the volume as the point included in all views
% orthoview(vol,centre) as above but allows the user to specify the centre
%  coordinate (i.e. centre must of the form [x y z] in voxels)
% orthoview(vol,centre,[imin imax]) as above but specifies the absolute 
%  intensity range to use imin <= intensity <= imax.   
%  This is useful for bypassing the default robust range when it is 
%  inappropriate.

if ((nargin>1) & (length(centre)>=3)),
  cvox=centre;
else
  cvox = round(size(vol)/2);
end

if ((nargin>2) & (length(voxdim)>=3)),
  vdims = abs(voxdim);
  vdims = vdims + (vdims==0)*1;  % get rid of zero dimensions
else
  vdims = [1 1 1];
end

if ((nargin>3) & (length(opts)>=2)),
  irange = [opts(1) opts(2)];
elseif ((nargin>3) & (length(opts)==1)),
  % Use absolute range
  irange = [min(min3(vol)) max(max3(vol))];
else
  % Use the robust min and max
  irange = percentile(vol,[2 98]);
end

if ( abs(irange(2) - irange(1))==0),
  % Use absolute range if robust is no good
  irange = [min(min3(vol)) max(max3(vol))];
end
if ( abs(irange(2) - irange(1))==0),
  irange = [irange(1) irange(1)+1];
end

ioff = irange(1);
iscale = (length(colormap)-1)/(irange(2) - irange(1));

clf

% Sagittal
subplot(2,2,1)
sagim = flipdim(permute(squeeze(vol(cvox(1),:,:)),[2 1]),1);
image(max(min((sagim-ioff)*iscale + 1,length(colormap)),1));
pbaspect([vdims(3) vdims(2) 1]);
axis off

% Coronal
subplot(2,2,2)
corim = flipdim(permute(squeeze(vol(:,cvox(2),:)),[2 1]),1);
image(max(min((corim-ioff)*iscale + 1,length(colormap)),1));
pbaspect([vdims(3) vdims(1) 1]);
axis off

% Axial
subplot(2,2,4)
axim = flipdim(permute(squeeze(vol(:,:,cvox(3))),[2 1]),1);
image(max(min((axim-ioff)*iscale + 1,length(colormap)),1));
pbaspect([vdims(2) vdims(1) 1]);
axis off

                                                                                                                                                 ./ov.m                                                                                              000664  002210  001747  00000000374 07330044540 012637  0                                                                                                    ustar 00mark                            analysis                        000000  000000                                                                                                                                                                         function ov(vol,centre,voxdims,opts)
% Calls orthoview

if (nargin>=4),
  orthoview(vol,centre,voxdims,opts);
elseif (nargin==3),
  orthoview(vol,centre,voxdims);
elseif (nargin==2),
  orthoview(vol,centre);
elseif (nargin==1),
  orthoview(vol);
end


                                                                                                                                                                                                                                                                    ./percentile.m                                                                                      000644  002206  001747  00000000470 07340176543 015253  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function [y] = percentile(x,p)
% y = percentile(x,p)
%   p is a number between 0 and 100 (the percentile value)
%   x is the data vector (it will reshape any matrix into a single column)

if (p>100), p=100; end
if (p<0),   p=0;   end
n = prod(size(x));
xx = sort(reshape(x,n,1));
y = xx(1 + round((n-1)*p/100));
                                                                                                                                                                                                        ./plot_mm.m                                                                                         000666  002206  001747  00000002356 07541622350 014574  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function plot_mm(mix,data,mm)

% plot_mm(mix,data)
% 
% plots Gaussian Mixture Model compared with data

nmixes = length(mix.centres);
res = (max(data) - min(data))/(length(data)/20);

mind = res*floor(min(data)/res);
maxd = res*ceil(max(data)/res);
al = [mind:res:maxd]';
alneg = [mind:res:-res]';
alplus = [res:res:maxd]';

mixtures = zeros(nmixes,length(al));

hs = hist(data,al);

bar(al,hs/sum(hs));

hold on;

for m=1:nmixes,
  if mix.priors(m) > 0,
    if(mm == 2 & m==3)      
      me=mix.centres(m);
      va=mix.covars(m);
      plus = mix.priors(m)*gammapdf(me^2/va, me/va, alplus)';
      mixtures(m,length(al)-length(alplus)+1:length(al)) = plus; 
      
    elseif(mm == 2 & m==1)
      me=-mix.centres(m);
      va=mix.covars(m);      
      neg = mix.priors(m)*gammapdf(me^2/va, me/va, -alneg)';
      mixtures(m,1:length(alneg)) = neg;      
      
    elseif(mix.covars(m)>0)
      mixtures(m,:) = mix.priors(m)*normpdf(al,mix.centres(m),sqrt(mix.covars(m)))';    
    end;
  end;
end;


for m=1:nmixes,
   
  %plot(al',mix.priors(m)*mixtures(m,:)./sum(sum(mixtures,2)),'r','LineWidth',2);
  plot(al',mixtures(m,:)./sum(sum(mixtures,2)),'r','LineWidth',2);
end;

plot(al',sum(mixtures,1)./sum(sum(mixtures,2)),'y--','LineWidth',2);

hold off;
                                                                                                                                                                                                                                                                                  ./plot_rms.m                                                                                        000664  002226  001747  00000000511 07323360061 014677  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function plot_rms(res)
% plot_rms(res)

clf
imno=1:1/6:(19-1/6);
hold on                   
axis([0 20 0 25])
for l=1:19,
  plot([l,l],[0 140],':');   
end
box on
set(gca,'FontSize',14);
xlabel('Image Number')
ylabel('RMS Deviation (mm)')
%title('FLIRT - Consistency Study')
h=gca;
set(h,'XTick',(1:18))
plot(imno,res(:,7),'-o')
                                                                                                                                                                                       ./plotseries.m                                                                                      000664  002226  001747  00000000162 07323360061 015233  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function plotseries(tseries, i, j, k, varargin)

plot(tseries);

title(sprintf('Slice %d; X: %d  Y: %d',k,i,j)); 
                                                                                                                                                                                                                                                                                                                                                                                                              ./range3.m                                                                                          000664  002210  001747  00000000115 07365553331 013376  0                                                                                                    ustar 00mark                            analysis                        000000  000000                                                                                                                                                                         function [r] = range3(x)
% returns [min3(x) max3(x)]

r = [min3(x) max3(x)];
                                                                                                                                                                                                                                                                                                                                                                                                                                                   ./read_avw.m                                                                                        000664  002226  001747  00000001345 07323360061 014636  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [img,dims,scales,bpp,endian] = read_avw(fname)
% [img, dims,scales,bpp,endian] = READ_AVW(fname)
%
%  Read in an analyse file into either a 3D or 4D
%  array (depending on the header information)
%  Ouput coordinates are in MEDx convention
%  except that all dimensions start at 1 rather than 0
%  Note: automatically detects char, short, long or double formats
%  Extracts the 4 dimensions (dims), 
%  4 scales (scales) and bytes per pixel (bpp) for voxels 
%  contained in the Analyse header file (fname)
%  Also returns endian = 'l' for little-endian or 'b' for big-endian
%
%  See also: READ_AVW_HDR, READ_AVW_IMG, SAVE_AVW, SAVE_AVW_HDR, SAVE_AVW_IMG

  img=read_avw_img(fname);
  [dims,scales,bpp,endian]= read_avw_hdr(fname);                                                                                                                                                                                                                                                                                           ./read_avw_complex.m                                                                                000664  002210  001747  00000002114 07325343475 015541  0                                                                                                    ustar 00mark                            analysis                        000000  000000                                                                                                                                                                         function [img,dims,scales,bpp,endian] = read_avw_complex(fname)
% [img, dims,scales,bpp,endian] = READ_AVW_COMPLEX(fname)
%
%  Read in a complex analyse file into either a 3D or 4D
%  array (depending on the header information)
%  Uses avwcomplex to make temporary files for reading
%  Ouput coordinates are in MEDx convention
%  except that all dimensions start at 1 rather than 0
%  Note: automatically detects char, short, long or double formats
%  Extracts the 4 dimensions (dims), 
%  4 scales (scales) and bytes per pixel (bpp) for voxels 
%  contained in the Analyse header file (fname)
%  Also returns endian = 'l' for little-endian or 'b' for big-endian
%
%  See also: READ_AVW_HDR, READ_AVW_IMG, SAVE_AVW, SAVE_AVW_HDR, SAVE_AVW_IMG

command=sprintf('! avwcomplex -realcartesian %s %s %s \n',fname,[fname,'R'],[fname,'I']);
eval(command);

[imgr,dims,scales,bpp,endian]=read_avw([fname,'R']);
[imgi,dims,scales,bpp,endian]=read_avw([fname,'I']);

img = imgr + j * imgi;

command=sprintf('! rm %s.hdr %s.img %s.hdr %s.img \n',[fname,'R'],[fname,'R'],[fname,'I'],[fname,'I']);
eval(command);
                                                                                                                                                                                                                                                                                                                                                                                                                                                    ./read_avw_hdr.m                                                                                    000664  002226  001747  00000002517 07323360061 015475  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [dims,scales,bpp,endian,datatype] = read_avw_hdr(fname)
% [dims,scales,bpp,endian] = READ_AVW_HDR(fname)
%
%  Extracts the 4 dimensions (dims), 
%   4 scales (scales) and bytes per pixel (bpp) for voxels 
%   contained in the Analyse header file (fname)
%   Also returns endian = 'l' for little-endian or 'b' for big-endian
%
%  See also: READ_AVW, READ_AVW_IMG, SAVE_AVW, SAVE_AVW_HDR, SAVE_AVW_IMG

% remove extension if it exists
if ( (length(findstr(fname,'.hdr'))>0) | ...
	(length(findstr(fname,'.img')>0)) ),
  fname=fname(1:(length(fname)-4));
end
fnhdr=strcat(fname,'.hdr');

% open file in big-endian
endian='b';
fid=fopen(fnhdr,'r','b');
testval = fread(fid,1,'int32');
% check if this gives the correct header size - if not use little-endian
if (testval~=348),
  fclose(fid);
  fid=fopen(fnhdr,'r','l');
  endian='l';
  testval = fread(fid,1,'int32');
  if (testval~=348),
    disp('Can not read this file format');
    return;
  end
end
	% ditch the remaining initial header stuff
  dummy=fread(fid,36,'char');
	% ditch dim[0] = No. dimensions
  dummy=fread(fid,1,'int16');
  dims=fread(fid,4,'int16');
  dummy=fread(fid,3,'int16');
  dummy=fread(fid,14,'char');
  datatype=fread(fid,1,'int16');  
  bpp=fread(fid,1,'int16');
  dummy=fread(fid,2,'char');
  dummy=fread(fid,1,'float');
  scales=fread(fid,4,'float');
fclose(fid);
return;
                                                                                                                                                                                 ./read_avw_img.m                                                                                    000664  002226  001747  00000002736 07323360061 015477  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function img = read_avw_img(filename);
%  [img] = READ_AVW_IMG(filename)
%
%  Read in an analyse file into either a 3D or 4D
%   array (depending on the header information)
%  Ouput coordinates are in MEDx convention
%   except that all dimensions start at 1 rather than 0
%  Note: automatically detects char, short, long or double formats
%  
%  See also: READ_AVW, READ_AVW_HDR, SAVE_AVW, SAVE_AVW_HDR, SAVE_AVW_IMG

% remove extension if it exists
if ( (length(findstr(filename,'.hdr'))>0) | ...
	(length(findstr(filename,'.img')>0)) ),
  filename=filename(1:(length(filename)-4));
end
fnimg=strcat(filename,'.img');
fnhdr=strcat(filename,'.hdr');

[dims,scales,bpp,endian,datatype] = read_avw_hdr(fnhdr);
fp=fopen(fnimg,'r',endian);
if (datatype==4),
  dat=fread(fp,'short');
elseif (datatype==2),
  dat=fread(fp,'char');
elseif (datatype==8),
  dat=fread(fp,'int');
elseif (datatype==64),
  dat=fread(fp,'double');
elseif (datatype==16),
   dat=fread(fp,'float32');
end
fclose(fp);

nvox = prod(dims);
if (length(dat)<nvox),
  error('Cannot open image as .img file does not contain as many voxels as the .hdr specifies');
elseif (length(dat)>nvox),
  disp('WARNING::truncating .img data as it contains more voxels than specified in the .hdr');
  dat = dat(1:nvox);
end


if (dims(4)>1),
  img = reshape(dat,dims(1),dims(2),dims(3),dims(4));
else
  img = reshape(dat,dims(1),dims(2),dims(3));
end

clear dat;

%% DEFUNCT FLIPPING
%% flip y dimension to be consistent with MEDx
%img=flipdim(img,2);
                                  ./read_dof.m                                                                                        000664  002226  001747  00000001202 07323360061 014601  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [dof] = read_dof(fname)
% [dof] = READ_DOF(fname)
%
%  Extracts the DOF list as a vector
%   from a .dof file generated by geomcal or mpr (both UMDS)
%

fid=fopen(fname);
while (~feof(fid)),
  str=fscanf(fid,'%s',1);
  if (strcmp(str,'DOF:')),
    dof=get_matrix(fid);
  end
end
fclose(fid);
return;


%-----------------------------------------

function g = get_matrix(fid)

num=fscanf(fid,'%d',1);
if (num<6)
 disp('Error in reading Linear_Transform matrix');
end
 g = fscanf(fid,'%f',inf);
 gl = length(g);
 g = reshape(g,gl/num,num).'; 
 gsz = size(g);
  % extract the last column
 if (gsz(2)>1)
   g = g(:,gsz(2));
 end
return;
                                                                                                                                                                                                                                                                                                                                                                                              ./read_medx_xfm.m                                                                                   000664  002226  001747  00000002515 07323360061 015650  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [g,a,s] = read_medx_xfm(fname)
% [g,a,s] = READ_MEDX_XFM(fname)
%
%  Extracts the GenericReslice matrix (g),
%   the AlignLinearReslice matrix (a), and
%   the ShadowTransform matrix (s)  (or UserTransform) from a
%   MEDx transformation file (given by fname)
%  Any 4x4 or 3x4 matrices should be formatted
%  Other, unrecognised, sizes will be an nx1 vector
%
g=[];
a=[];
s=[];
fid=fopen(fname);
if (fid<0),
  error(['File ',fname,' not found.']);
end
while (~feof(fid)),
  str=fscanf(fid,'%s',1);
  if (strcmp(str,'/GenericReslice')),
    tstg=get_matrix(fid);
    if length(tstg>0),  g=tstg;  end
  end
  if (strcmp(str,'/MotionCorrectionReslice')),
    tsta=get_matrix(fid);
    if length(tsta>0),  a=tsta;  end
  end
  if (strcmp(str,'/ShadowTransform')),
    tsts=get_matrix(fid);
    if length(tsts>0),  s=tsts;  end
  end
end
fclose(fid);
return;


%-----------------------------------------

function g = get_matrix(fid)

g=[];
str1 = fscanf(fid,'%s',2);
if (~strcmp(str1,'<</matrix')),
  return;
end
str1 = fscanf(fid,'%s',1);
if (~strcmp(str1,'[')),
 disp('Error in reading Reslice Matrix');
end
 str1 = fscanf(fid,'%s',1);
 while (~strcmp(str1,']')),
   g=[g str2num(str1)];
%   disp(['READ ',str1]);
   str1 = fscanf(fid,'%s',1);
 end
if (length(g)==16),
  g=reshape(g,4,4).';
elseif (length(g)==12),
  g=reshape(g,4,3).';
end
return;
                                                                                                                                                                                   ./read_mri_xfm.m                                                                                    000664  002226  001747  00000001160 07323360061 015475  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [transmat] = read_mri_xfm(fname)
% [transmat] = READ_MRI_XFM(fname)
%
%  Extracts the 4x4 transformation matrix
%   from a .xfm file generated by minctracc (via mritotal)
%

fid=fopen(fname);
while (~feof(fid)),
  str=fscanf(fid,'%s',1);
  if (strcmp(str,'Linear_Transform')),
    transmat=get_matrix(fid);
  end
end
fclose(fid);
return;


%-----------------------------------------

function g = get_matrix(fid)

str1 = fscanf(fid,'%s',1);
if (~strcmp(str1,'=')),
 disp('Error in reading Linear_Transform matrix');
end
 g = fscanf(fid,'%f',12);
 g = reshape(g,1,12);
 g = [g 0 0 0 1];
 g=reshape(g,4,4).';
return;
                                                                                                                                                                                                                                                                                                                                                                                                                ./read_vest.m                                                                                       000644  002206  001747  00000001167 07327517741 015105  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function x = read_vest(filename)

fp=fopen(filename,'r');

while(1)
	line = fgetl(fp);
	if(strcmp(line,'/Matrix'))
		break;
	end;
end;

str = fgetl(fp);
x = read_vest_line(str);
str = fgetl(fp);

c=1;
while(isstr(str))
	c=c+1;
	x(c,:) = read_vest_line(str);
	str = fgetl(fp);
end;

fclose(fp);

%%%%%%%%%%%%%%%%

function val = read_vest_line(str)

indDelim = sort([0,findstr(str,sprintf('\t')),findstr(str,sprintf(' '))]);
val = zeros(1,length(indDelim)-1);
s=0;

for ctDelim = 1:length(indDelim)-1,
      s=s+1;
      val(s) = str2double(str(indDelim(ctDelim)+1:indDelim(ctDelim+1)-1));
end

if any(isnan(val)),
   val=NaN;
end

                                                                                                                                                                                                                                                                                                                                                                                                         ./save_avw.m                                                                                        000664  002226  001747  00000001235 07323360061 014657  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function save_avw(img,fname,vtype,vsize)
% SAVE_AVW(img,fname,vtype,vsize) 
%
%  Create and save an analyse header (.hdr) and image (.img) file
%   for either a 2D or 3D or 4D array (automatically determined).
%  Note: 
%        assumes the data uses a MEDx coordinate convention
%	 and that slice direction is transverse
%   
%  vtype is a single character string: 'b' byte, 's' short or 'f' float
%  vsize is a vector [x y z tr] containing the voxel sizes in mm and
%  the tr in seconds  (defaults: [1 1 1 3])
%
%  See also: SAVE_AVW_HDR, SAVE_AVW_IMG, READ_AVW, READ_AVW_HDR, READ_AVW_IMG
%

   save_avw_hdr(img,fname,vtype,vsize);
   save_avw_img(img,fname,vtype);
                                                                                                                                                                                                                                                                                                                                                                   ./save_avw_complex.m                                                                                000664  002210  001747  00000001677 07325342675 015602  0                                                                                                    ustar 00mark                            analysis                        000000  000000                                                                                                                                                                         function save_avw_complex(img,fname,vsize)
% SAVE_AVW_COMPLEX(img,fname,vsize) 
%
%  Create and save an analyse header (.hdr) and image (.img) file
%   for either a 2D or 3D or 4D array (automatically determined).
%  Only for use with complex data.  (uses avwcomplex)
%
%  Note: 
%        assumes the data uses a MEDx coordinate convention
%	 and that slice direction is transverse
%   
%  vtype is a single character string: 'b' byte, 's' short or 'f' float
%  vsize is a vector [x y z tr] containing the voxel sizes in mm and
%  the tr in seconds  (defaults: [1 1 1 3])
%
%  See also: SAVE_AVW_HDR, SAVE_AVW_IMG, READ_AVW, READ_AVW_HDR, READ_AVW_IMG
%

save_avw(real(img),[fname,'R'],'f',vsize);
save_avw(imag(img),[fname,'I'],'f',vsize);
command=sprintf('! avwcomplex -complex %s %s %s \n',[fname,'R'],[fname,'I'],fname);
eval(command);
command=sprintf('! rm %s.hdr %s.img %s.hdr %s.img \n',[fname,'R'],[fname,'R'],[fname,'I'],[fname,'I']);
eval(command);
                                                                 ./save_avw_hdr.m                                                                                    000664  002226  001747  00000004352 07326060214 015517  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function save_avw_hdr(img,fname,vtype,vsize)
% SAVE_AVW_HDR(img,fname,vtype,vsize) 
%
%  Create and save an analyse header file
%   for either a 2D or 3D or 4D array (automatically determined).
%  Note: 
%        assumes the data uses a MEDx coordinate convention
%	 and that slice direction is transverse
%   
%  vtype is a single character string: 'b' byte, 's' short or 'f' float
%  vsize is a vector [x y z tr] containing the voxel sizes in mm and
%  the tr in seconds  (defaults: [1 1 1 3])
%
%  See also: SAVE_AVW, SAVE_AVW_IMG, READ_AVW, READ_AVW_HDR, READ_AVW_IMG 
%

% swap first and second argument in case save_avw_img convention is
% used
check=length(size(fname));
if(check~=2)
   tmp=img;
   img=fname;
   fname=tmp;
end

% remove extension if it exists
if ( (length(findstr(fname,'.hdr'))>0) | ...
        (length(findstr(fname,'.img')>0)) ),
  fname=fname(1:(length(fname)-4));
end

% remove input file:
!touch hdrcreate.txt;
!rm hdrcreate.txt;

% remove headerfile
fname2=strcat(fname,'.hdr');
tmpstr1=sprintf('!touch %s', fname2);
tmpstr2=sprintf('!rm %s', fname2);

eval(tmpstr1);
eval(tmpstr2);

% establish dynamic range
imgmax=ceil(max(max(max(max(img)))));
imgmin=floor(min(min(min(min(img)))));

% create file to use as input into header program
dims = [size(img) 1 1];

if(nargin==2)
  vtype='s';
  vsize=[1 1 1 3];
elseif(nargin==3)
  tmp=size(vtype);
  if(tmp(2)==1)
     vsize=[1 1 1 3];
  else
     vsize=vtype;
     if size(vsize,2)==3
	vsize=[vsize 3];
     end;
     vtype='s';
  end
else
  tmp=size(vtype);
  if(tmp(2)==3)
     tmp2=vtype;
     vtype=vsize;
     vsize=tmp2;
  end
end

if (length(vsize)<3),
  vsize(3)=1;
end
if (length(vsize)<4),
  vsize(4)=3;
end

 fid = fopen('hdrcreate.txt','w');
  if imgmin~=imgmax
    fprintf(fid,'%d\n%d\n%d\n%d\n%s\nv\n%6.4f\n%6.4f\n%6.4f\n%6.4f\nr\n%6.0f\n%6.0f\ns\n',dims(1),dims(2), dims(3),dims(4),vtype,vsize(1),vsize(2),vsize(3),vsize(4),imgmin,imgmax);
  else
    fprintf(fid,'%d\n%d\n%d\n%d\n%s\nv\n%6.4f\n%6.4f\n%6.4f\n %6.4f\ns\n',dims(1),dims(2), dims(3),dims(4),vtype,vsize(1),vsize(2),vsize(3), vsize(4));  
  end;
  fclose(fid);

% call header program

tmp=sprintf('! /usr/local/fsl/bin/header -n %s < hdrcreate.txt \n',fname);
eval(tmp);
disp(' ');

% remove input file:
!rm hdrcreate.txt;

                                                                                                                                                                                                                                                                                      ./save_avw_img.m                                                                                    000664  002226  001747  00000002131 07323360061 015507  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function save_avw_img(img,fname,vtype);
%  SAVE_AVW_IMG(img,fname,vtype)
%
%  Save an array (img) as an analyse file (only the .img) 
%   for either a 2D or 3D or 4D array (automatically determined)
%  Note: 
%        assumes the data uses a MEDx coordinate convention
% 
%  vtype is a single character string: 'b' byte, 's' short or 'f' float
%
%  See also: SAVE_AVW, SAVE_AVW_HDR, READ_AVW, READ_AVW_HDR, READ_AVW_IMG
%

% swap first and second argument in case save_avw_img convention is
% used
check=length(size(fname));
if(check~=2)
   tmp=img;
   img=fname;
   fname=tmp;
end

% remove extension if it exists
if ( (length(findstr(fname,'.hdr'))>0) | ...
	(length(findstr(fname,'.img')>0)) ),
  fname=fname(1:(length(fname)-4));
end
fnimg=strcat(fname,'.img');

fp=fopen(fnimg,'w');
dims = size(img);

%% DEFUNCT
%% flip y dimension to be consistent with MEDx
%% dat=flipdim(img,2);

dat = img;
dat = reshape(dat,prod(dims),1);

switch vtype
  case 'f'
    vtype2='float';
  case 's'
    vtype2='short';
  case 'b'
    vtype2='schar';
  case 'u'
    vtype2='uchar';
end;

fwrite(fp,dat,vtype2);
fclose(fp);

                                                                                                                                                                                                                                                                                                                                                                                                                                       ./save_vest.m                                                                                       000644  002206  001747  00000000712 07340202406 015103  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function write_vest(x,filename)

% write_vest(x,filename)
%
% writes x in vest format where size(x) is ntpts*nevs (NumPoints*NumWaves)

fp=fopen(filename,'w');

dims = size(x);

fprintf(fp,'! VEST-Waveform File\n');
fprintf(fp,'/NumWaves\t%i\n',dims(2));
fprintf(fp,'/NumPoints\t%i\n',dims(1));
fprintf(fp,'/Skip\n');
fprintf(fp,'\n/Matrix\n');

for t=1:dims(1),
  for e=1:dims(2),
    fprintf(fp,'%e\t',x(t,e));
  end;
  fprintf(fp,'\n');
end;

fclose(fp);
                                                      ./scaleim.m                                                                                         000664  002226  001747  00000000557 07327575324 014505  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function imsc = scaleim(im)
% imsc = scaleim(im)
% 
% Returns a scaled version of the image, with intensities between 0 and 255

imsc = min(255,(im - min(min(im)))/(max(max(im)) - min(min(im))) * 256);

% Use the robust min and max
irange = percentile(im,[2 98]);
ioff = irange(1);
iscale = 255/(irange(2) - irange(1) + 1);

imsc = clamp((im - ioff)*iscale, 0,255);

                                                                                                                                                 ./sigmoid.m                                                                                         000644  017223  001747  00000000052 07363004301 014340  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function b=sigmoid(a);
b=1./(1+exp(-a));

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ./spm2medx.m                                                                                        000664  002226  001747  00000001056 07323360061 014604  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [transmat] = spm2medx(tr_mat,init_size,final_size)
% [transmat] = spm2medx(tr_mat,init_size,final_size)
% Converts a transformation matrix between spm format
%   and medx format  (pre and post flip in x and y)
% It must be given the sizes of the initial and final
%   volumes init_size, final_size both as [xmax,ymax,zmax].
%   Convention is that xmax = no x voxels , ymax = etc.

flip1=[-1 0 0 init_size(1); 0 -1 0 init_size(2); 0 0 1 1; 0 0 0 1];
flip2=[-1 0 0 final_size(1); 0 -1 0 final_size(2); 0 0 1 -1; 0 0 0 1];
transmat=flip2*tr_mat*flip1;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ./squash.m                                                                                          000664  002206  001747  00000000265 07403713217 014424  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function [ret, dims] = squash(x)

% function [ret, dims] = squash(x)
%
% performs ret = reshape(x,prod(size(x)),1);
% see unsquash

dims = size(x);
ret = reshape(x,prod(size(x)),1);                                                                                                                                                                                                                                                                                                                                           ./stochastic_design.m                                                                               000644  002206  001747  00000003177 07521765400 016622  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function [d, t] = stochastic_design(tr, meanisi, minisi, res, n, isidist, threecolfile) 

% function [d, t] = stochastic_design(tr, meanisi, minisi,
% resolution, numscans, isidist, threecolfile)
%
% tr, meanisi, minisi, resolution in time units.
%
% n is the number of scans
%
% isidist should be 3 for fixed isi designs
%
% isidist should be 2 for poisson isi designs
%
% isidist should be 1 for uniform isi designs, the
% isi's are sampled from uniform distribution U(meanisi-sqrt(meanisi), meanisi+sqrt(meanisi))
% and no values less then minisi get used.
%
% isidist should be 0 for normal distribution isi designs
% the isi's are sampled from N(meanisi, meanisi)
% and no values less then minisi get used.
%
% Returns vector d of ones and zeros and time vector t
% and efficiency = 1/variance(b) = 1/inv(d'*d);
% If threecolfile ~= '' then writes out the design in FEAT's 3 column format
% to a file specified by threecolfile.

if isidist == 0,
  rs = randn(round(4*n*tr/meanisi),1)*(sqrt(meanisi)) + meanisi;
elseif isidist == 1,
  rs = (rand(round(4*n*tr/meanisi),1)-0.5)*(sqrt(meanisi)) + meanisi;
elseif isidist == 2,
  rs = poissrnd(meanisi,round(4*n*tr/meanisi),1);
elseif isidist == 3,
  rs = ones(round(4*n*tr/meanisi),1)*meanisi;
end;

rs2 = rs(find(rs >= minisi));

crs2=cumsum(rs2);

[i,j]=find(crs2>tr*n);
%rs2(i(1))=tr*n-crs2(i(1)-1);
rs3=rs2(1:i(1)-1);

if(~strcmp(threecolfile,'')),
x=zeros(length(rs3),3);
x(:,1) = cumsum(rs3);
x(:,2)=0.1;
x(:,3)=1;

fid = fopen(threecolfile,'w');
fprintf(fid,'%6.4f  %6.2f %6.2f\n',x');
fclose(fid);
end;

t = [res:res:n*tr];

%keyboard;

d = zeros(n*tr/res,1);
crs3 = cumsum(rs3);
d(round(crs3/res)) = 1;
                                                                                                                                                                                                                                                                                                                                                                                                 ./subplotposition.m                                                                                 000664  002226  001747  00000001355 07323360061 016324  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [left, bottom, width, height] = subplotposition(x,numx,y,numy,spacex,spacey)

% [left, bottom, width, height] = subplotposition(x,numx,y,numy,spacex,spacey)
% 
% Calls subplot passing a calculated position
% numx is the number of columns, x is the column to use counted
% from the left
%
% numy is the number of rows, y is the row to use counted from the
% bottom
%
% spacex and spacey are optional and specify the spacing between plots 
% the default is 0.05

if nargin < 5,
  spacex = 0.05;
end;

if nargin < 6,
  spacey = spacex;
end;

height = (1-spacey)/numy - spacey;
bottom = (y-1)*(height+spacey) + spacey;

width = (1-spacex)/numx - spacex;
left = (x-1)*(width+spacex) + spacex;

subplot('position', [left, bottom, width, height]);                                                                                                                                                                                                                                                                                   ./subsample.m                                                                                       000664  002226  001747  00000000536 07323360061 015042  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [subim] = subsample(im)
% subim = subsample(im)
%
% Subsamples the image (im) by a factor of 2, with pre-blurring, and
%  returns the image as subim

[M,N] = size(im);
M2 = floor(M/2);
N2 = floor(N/2);
imblur = conv2(im,gauss(0.85,3),'same');
subim = zeros(M2,N2);
for x1=1:M2,
 for y1=1:N2,
  subim(x1,y1) = imblur(x1*2-1,y1*2-1);
 end
end
                                                                                                                                                                  ./sum3.m                                                                                            000664  002226  001747  00000000111 07323360061 013723  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [m] = sum3(x)
% performs sum(sum(sum(x)))

m = sum(sum(sum(x)));                                                                                                                                                                                                                                                                                                                                                                                                                                                       ./textscatter.m                                                                                     000644  017335  001765  00000000650 07513257264 014052  0                                                                                                    ustar 00cedric                          fsl                             000000  000000                                                                                                                                                                         function textscatter(a,b,style);
% TEXTSCATTER(A,B,STYLE)
%
% puts the index numbers next to the points in a scatter plot
% A,B are the 2 vectors
% Style is the third argument to plot - e.g 'ro'
%
% (C) Cedric - 2002

if(length(a)~=length(b))
     error('bugger off')
end
if(nargin<3)
     style='ro';
end
hold on
for i=1:length(a)
     h=plot(a(i),b(i),style);
     A=get(h);
     text(A.XData+0.03,A.YData,num2str(i))
end
                                                                                        ./tfitseries.m                                                                                      000644  002206  001747  00000000231 07450622762 015276  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function tfitseries(tseries, i, j, k, varargin)

if std(tseries)<=0,return;end;

multitfit(tseries,1);

title(sprintf('Slice %d; X: %d  Y: %d',k,i,j)); 
                                                                                                                                                                                                                                                                                                                                                                       ./threecol2ev.m                                                                                     000666  002203  001747  00000001513 07515560112 014570  0                                                                                                    ustar 00heidi                           analysis                        000000  000000                                                                                                                                                                         function ev = threecol2ev(t,res,nsecs);
%THREECOL2EV convert 3col format to ev format.
% EV = THREECOL2EV(T,RES);
%
% RES is required temporal resolution
%
% T is the 3 column format
% 
% EV = THREECOL2EV(T,RES,NSECS);
%
% NSECS is the total length in seconds of your data.  
% If this argument is not included then total length 
% is set to the length of your three column ev
% 
% TB, HJB
if(nargin<2)
  error('Not enough input arguments')
end
t_res=t;
t_res(:,1:2)=round(t(:,1:2)/res);
if(nargin==2) 
  ev = zeros(t_res(size(t_res,1),1)+t_res(size(t_res,1),2),1);
elseif(nargin==3)
  if(nsecs<(t(size(t,1),1)+t(size(t,1),2)))
    error('Your three column ev is longer than nsecs')
  else
    ev = zeros(nsecs/res,1);
  end
end
%t_res(:,1)=t_res(:,1)+1-t_res(1,1);
for i=1:size(t_res,1)
  ev(t_res(i,1):t_res(i,1)+t_res(i,2))=t_res(i,3);
end
                                                                                                                                                                                     ./thresh.m                                                                                          000644  002206  001747  00000000216 07351424602 014406  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function x = thresh(x, t)

% x = thresh(x, t)
% set all x>t to t and all x<t to 0

x2 = squash(x);
x2(find(x<t)) = 0;
x = reshape(x2,size(x));                                                                                                                                                                                                                                                                                                                                                                                  ./unsquash.m                                                                                        000644  002206  001747  00000000212 07403713363 014757  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function [ret] = unsquash(x,dims)

% function [ret] = unsquash(x)
%
% performs ret = reshape(x,dims);
% see squash

ret = reshape(x,dims);                                                                                                                                                                                                                                                                                                                                                                                      ./var3.m                                                                                            000664  002226  001747  00000000133 07323360061 013713  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function [v] = var3(x)
% performs mean3((x - mean3(x)).^2);

v = mean3((x - mean3(x)).^2);
                                                                                                                                                                                                                                                                                                                                                                                                                                     ./vbmix2mix.m                                                                                       000644  002206  001747  00000000503 07421036015 015030  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function mix = vbmix2mix(vbmix)

% mix = vbmix2mix(vbmix)
%
% converts Tim's dodgey Gaussian mixture struct to the
% proper Netlab struct.

mix = gmm(1,vbmix.ncentres,'spherical');
mix.priors = vbmix.mixers.lambda/sum(vbmix.mixers.lambda);
mix.centres = vbmix.centres.means';
mix.covars = 1./(vbmix.precs.b.*vbmix.precs.c);                                                                                                                                                                                             ./view3d.m                                                                                          000644  017223  001747  00000001710 07345457335 014133  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function view3d(vol,x,y,z);
% VIEW3D - Show 3D sliced view of a volume
%
%  VIEW3D(VOL,X,Y,Z)
%  VOL is volume
%  X is sagittal slice number
%  Y is coronal slice number
%  Z is axial slice number
%
% copyright (C) T.Behrens 2001

sagittal=squeeze(vol(x,:,:));
coronal=squeeze(vol(:,y,:));
axial=squeeze(vol(:,:,z));

sag=surf(sagittal/1000000+x,sagittal,'linestyle','none');
FVCsag=surf2patch(sag);
FVCsag1=FVCsag;
FVCsag1.vertices(:,1)=FVCsag.vertices(:,3);
FVCsag1.vertices(:,3)=FVCsag.vertices(:,1);
delete(sag);

cor=surf(coronal'/1000000+y,coronal','linestyle','none');
FVCcor=surf2patch(cor);
FVCcor1=FVCcor;
FVCcor1.vertices(:,2)=FVCcor.vertices(:,3);
FVCcor1.vertices(:,3)=FVCcor.vertices(:,2);
delete(cor);

ax=surf(axial'/1000000+z,axial','linestyle','none');
FVCax=surf2patch(ax);
delete(ax);

colormap('bone'); hold on;
patch(FVCsag1,'linestyle','none');
patch(FVCcor1,'linestyle','none');
patch(FVCax,'linestyle','none');
shading faceted; view(3)







                                                        ./view_feat_results.m                                                                               000664  002206  001747  00000002715 07522267227 016663  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function view_feat_results(featdir, copenum, t_or_f, zimg)

% view_feat_results(featdir, c)
%
% featdir is the feat directory.
% runs ViewFMRI using glm_fit 
% to show univariate analysis
% z-statistics and model fit for
% t-contrast number c.
%
% view_feat_results(featdir, c, t_or_f)
%
% shows model fit for the:
% t-contrast number c if t_or_f = 't'
% f-contrast number c if t_or_f = 'f'
%
% view_feat_results(featdir, c, t_or_f, zimg) 
%
% if zimg = 0 shows clustered and rendered zstats
% if zimg = 1 shows univariate analysis zstats

global xs dm copes zs acs c fc;

if (nargin <4) zimg=1; end;
if (nargin <3) t_or_f='t'; end;

statsdir = strcat(featdir,'/stats');
xs = read_avw(sprintf('%s/filtered_func_data', featdir));

if zimg==1,
  if(strcmp(t_or_f,'f'))
    zs = read_avw(sprintf('%s/zfstat%d', statsdir,copenum));
  else
    zs = read_avw(sprintf('%s/zstat%d', statsdir,copenum));
  end;
else,
  if(strcmp(t_or_f,'f'))
    zs = read_avw(sprintf('%s/Rendered_thresh_zfstat%d', ...
			      featdir,copenum));
  else
    zs = read_avw(sprintf('%s/Rendered_thresh_zstat%d', ...
			  featdir,copenum));
  end;
end;

copes = read_avw(sprintf('%s/cope%d', statsdir,copenum));
dm = read_vest(sprintf('%s/design.mat', featdir));

if(strcmp(t_or_f,'t'))
  c = read_vest(sprintf('%s/design.con', featdir));
else,
  fc = read_vest(sprintf('%s/design.fts', featdir));
end;

acs = read_avw(sprintf('%s/threshac1',statsdir));

ViewFMRI(xs, zs, 'glm_fit', 'glm_fit', copenum, t_or_f);

                                                   ./when.m                                                                                            000644  017223  001747  00000001042 07374446455 013674  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function when;
fid=fopen('thisisme');
if(fid == -1)
  !whoami > thisisme;
  fid=fopen('thisisme');
  S=fscanf(fid,'%s');
  !rm thisisme
  if(strcmp(S,'beckmann'))
    disp('1945');
  elseif(strcmp(S,'woolrich'))
    disp('2003 if you are lucky');
  elseif(strcmp(S,'mark'))
    disp('Why does it matter? You will be in Australia');
  elseif(strcmp(S,'prb'))
    disp('You were rowing mate');  
  elseif(strcmp(S,'steve'))
    disp('FSL indeed - You still need MATLAB. May as well use SPM');   
  else
    disp(date)
  end
else
  disp(date)
end


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              ./where.m                                                                                           000644  017223  001747  00000000047 07352653237 014043  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function where;
disp('Under your nose')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         ./write_medx_xfm.m                                                                                  000664  002226  001747  00000004200 07323360061 016060  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function write_medx_xfm(fname,no_voxels,g,a,s,outputdim)
% WRITE_MEDX_XFM(fname,no_voxels,g,a,s,outputdim)
%
%  Writes a MEDx transformation file using
%   the GenericReslice matrix (g),
%   the AlignLinearReslice matrix (a), 
%   the ShadowTransform matrix (s).
%   into a file specified by fname.
%  Any 4x4 or 3x4 matrices will be correctly formatted
%  The outputdim matrix (4x4) gives the sizes for the voxel dimensions
%   e.g. outputdim = diag([0.93 0.93 5 1])
%  The no_voxels gives the volume size in voxels [no_x,no_y,no_z]
%   e.g. no_voxels = [256 256 30]
%
fid=fopen(fname,'w');
fprintf(fid,'%s\n','%!VEST-Transformations');
fprintf(fid,'%s\n','<<');
if (length(g)>0),
     fprintf(fid,'%s\n','    /GenericReslice     <<');
     put_matrix(fid,'matrix',g);
     fprintf(fid,'%s\n','        /order  1');
     put_matrix(fid,'outputsize',no_voxels);
     if (length(outputdim)>0),
      put_matrix(fid,'outputusermatrix',outputdim);
     end
     fprintf(fid,'%s\n','    >>');
end
if (length(a)>0),
     fprintf(fid,'%s\n','    /MotionCorrectionReslice     <<');
     put_matrix(fid,'matrix',a);
     fprintf(fid,'%s\n','        /order  1');
     put_matrix(fid,'outputsize',no_voxels);
     if (length(outputdim)>0),
      put_matrix(fid,'outputusermatrix',outputdim);
     end
     fprintf(fid,'%s\n','        /outputunits    (mm)');
     fprintf(fid,'%s\n','        /talairachcalibrated?   false');
     fprintf(fid,'%s\n','    >>');
end
if (length(s)>0),
     fprintf(fid,'%s\n','    /ShadowTransform     <<');
     put_matrix(fid,'matrix',s);
     fprintf(fid,'%s\n','        /order  1');
     put_matrix(fid,'outputsize',no_voxels);
     if (length(outputdim)>0),
      put_matrix(fid,'outputusermatrix',outputdim);
     end
     fprintf(fid,'%s\n','    >>');
end
fprintf(fid,'%s\n','>>');
fclose(fid);
return;


%-----------------------------------------

function put_matrix(fid,name,g)

if (length(g)<=0),
  return;
end
if (length(g)~=prod(size(g))),
  g=reshape(g.',prod(size(g)),1);
end
fprintf(fid,'%s\n',['        /',name,' [']);
for idx=1:length(g),
  fprintf(fid,'%s\n',['            ',num2str(g(idx))]);
end
fprintf(fid,'%s\n','        ]');
return;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ./view3d.m                                                                                          000644  017223  001747  00000001710 07345457335 014133  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function view3d(vol,x,y,z);
% VIEW3D - Show 3D sliced view of a volume
%
%  VIEW3D(VOL,X,Y,Z)
%  VOL is volume
%  X is sagittal slice number
%  Y is coronal slice number
%  Z is axial slice number
%
% copyright (C) T.Behrens 2001

sagittal=squeeze(vol(x,:,:));
coronal=squeeze(vol(:,y,:));
axial=squeeze(vol(:,:,z));

sag=surf(sagittal/1000000+x,sagittal,'linestyle','none');
FVCsag=surf2patch(sag);
FVCsag1=FVCsag;
FVCsag1.vertices(:,1)=FVCsag.vertices(:,3);
FVCsag1.vertices(:,3)=FVCsag.vertices(:,1);
delete(sag);

cor=surf(coronal'/1000000+y,coronal','linestyle','none');
FVCcor=surf2patch(cor);
FVCcor1=FVCcor;
FVCcor1.vertices(:,2)=FVCcor.vertices(:,3);
FVCcor1.vertices(:,3)=FVCcor.vertices(:,2);
delete(cor);

ax=surf(axial'/1000000+z,axial','linestyle','none');
FVCax=surf2patch(ax);
delete(ax);

colormap('bone'); hold on;
patch(FVCsag1,'linestyle','none');
patch(FVCcor1,'linestyle','none');
patch(FVCax,'linestyle','none');
shading faceted; view(3)







                                                        ./view_feat_results.m                                                                               000664  002206  001747  00000002715 07522267227 016663  0                                                                                                    ustar 00woolrich                        analysis                        000000  000000                                                                                                                                                                         function view_feat_results(featdir, copenum, t_or_f, zimg)

% view_feat_results(featdir, c)
%
% featdir is the feat directory.
% runs ViewFMRI using glm_fit 
% to show univariate analysis
% z-statistics and model fit for
% t-contrast number c.
%
% view_feat_results(featdir, c, t_or_f)
%
% shows model fit for the:
% t-contrast number c if t_or_f = 't'
% f-contrast number c if t_or_f = 'f'
%
% view_feat_results(featdir, c, t_or_f, zimg) 
%
% if zimg = 0 shows clustered and rendered zstats
% if zimg = 1 shows univariate analysis zstats

global xs dm copes zs acs c fc;

if (nargin <4) zimg=1; end;
if (nargin <3) t_or_f='t'; end;

statsdir = strcat(featdir,'/stats');
xs = read_avw(sprintf('%s/filtered_func_data', featdir));

if zimg==1,
  if(strcmp(t_or_f,'f'))
    zs = read_avw(sprintf('%s/zfstat%d', statsdir,copenum));
  else
    zs = read_avw(sprintf('%s/zstat%d', statsdir,copenum));
  end;
else,
  if(strcmp(t_or_f,'f'))
    zs = read_avw(sprintf('%s/Rendered_thresh_zfstat%d', ...
			      featdir,copenum));
  else
    zs = read_avw(sprintf('%s/Rendered_thresh_zstat%d', ...
			  featdir,copenum));
  end;
end;

copes = read_avw(sprintf('%s/cope%d', statsdir,copenum));
dm = read_vest(sprintf('%s/design.mat', featdir));

if(strcmp(t_or_f,'t'))
  c = read_vest(sprintf('%s/design.con', featdir));
else,
  fc = read_vest(sprintf('%s/design.fts', featdir));
end;

acs = read_avw(sprintf('%s/threshac1',statsdir));

ViewFMRI(xs, zs, 'glm_fit', 'glm_fit', copenum, t_or_f);

                                                   ./when.m                                                                                            000644  017223  001747  00000001042 07374446455 013674  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function when;
fid=fopen('thisisme');
if(fid == -1)
  !whoami > thisisme;
  fid=fopen('thisisme');
  S=fscanf(fid,'%s');
  !rm thisisme
  if(strcmp(S,'beckmann'))
    disp('1945');
  elseif(strcmp(S,'woolrich'))
    disp('2003 if you are lucky');
  elseif(strcmp(S,'mark'))
    disp('Why does it matter? You will be in Australia');
  elseif(strcmp(S,'prb'))
    disp('You were rowing mate');  
  elseif(strcmp(S,'steve'))
    disp('FSL indeed - You still need MATLAB. May as well use SPM');   
  else
    disp(date)
  end
else
  disp(date)
end


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              ./where.m                                                                                           000644  017223  001747  00000000047 07352653237 014043  0                                                                                                    ustar 00behrens                         analysis                        000000  000000                                                                                                                                                                         function where;
disp('Under your nose')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         ./write_medx_xfm.m                                                                                  000664  002226  001747  00000004200 07323360061 016060  0                                                                                                    ustar 00beckmann                        analysis                        000000  000000                                                                                                                                                                         function write_medx_xfm(fname,no_voxels,g,a,s,outputdim)
% WRITE_MEDX_XFM(fname,no_voxels,g,a,s,outputdim)
%
%  Writes a MEDx transformation file using
%   the GenericReslice matrix (g),
%   the AlignLinearReslice matrix (a), 
%   the ShadowTransform matrix (s).
%   into a file specified by fname.
%  Any 4x4 or 3x4 matrices will be correctly formatted
%  The outputdim matrix (4x4) gives the sizes for the voxel dimensions
%   e.g. outputdim = diag([0.93 0.93 5 1])
%  The no_voxels gives the volume size in vo