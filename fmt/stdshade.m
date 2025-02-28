function stdshade(amatrix,alpha,acolor,F,lw,smth,varargin)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end

if exist('F','var')==0 || isempty(F); 
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1;0;
end  

if ne(size(F,1),1)
    F=F';
end

if ~exist('lw'); lw=1; end

if length(varargin) == 0
    ls = '-'; 
else
    ls = varargin{1};  
end

amean=smooth(nanmean(amatrix),smth)';
% astd=nanstd(amatrix); % to get std shading
 astd=nanstd(amatrix)/sqrt(size(amatrix,1)); % to get sem shading
% astd = nanstd(amatrix) ./ sqrt(sum(~isnan(amatrix)));% to get sem shading iwht proper non nan value


if exist('alpha','var')==0 || isempty(alpha) 
    fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none');
    acolor='k';
else fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');    
end

if ishold==0
    check=true; else check=false;
end

hold on;plot(F,amean,'Color',acolor,'linewidth',lw,'linestyle',ls); %% change color or linewidth to adjust mean line

if check
    hold off;
end


end



