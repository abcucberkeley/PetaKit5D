function [window] = tukwin(data,r,useGPU)
    [sy,sx,sz]=size(data);
    
    if useGPU
        window=gpuArray(ones(sy,sx,sz));
    else
        window=ones(sy,sx,sz);
    end
    
    [xx,yy,zz]=meshgrid(0:1/(sx-1):1,0:1/(sy-1):1,0:1/(sz-1):1);
    
    
    index1x = (0<=xx)&(xx<=r/2);
    index2x = (r/2<=xx)&(xx<=(1-r/2));
    index3x = ((1-r/2)<=xx)&(xx<=1);
    
    index1y = (0<=yy)&(yy<=r/2);
    index2y = (r/2<=yy)&(yy<=(1-r/2));
    index3y = ((1-r/2)<=yy)&(yy<=1);
    
    if(useGPU)
        index1x = gpuArray(index1x);
        index2x = gpuArray(index2x);
        index3x = gpuArray(index3x);

        index1y = gpuArray(index1y);
        index2y = gpuArray(index2y);
        index3y = gpuArray(index3y);
    end
    
%     index1z = (0<=zz)&(zz<=r/2);
%     index2z = (r/2<=zz)&(zz<=(1-r/2));
%     index3z = ((1-r/2)<=zz)&(zz<=1);
    
    window(index1x)=min(window(index1x),0.5*(1+cos(2*pi/r*(xx(index1x)-r/2))));
    window(index2x)=min(window(index2x),1);
    window(index3x)=min(window(index3x),0.5*(1+cos(2*pi/r*(xx(index3x)-1+r/2))));

    window(index1y)=min(window(index1y),0.5*(1+cos(2*pi/r*(yy(index1y)-r/2))));
    window(index2y)=min(window(index2y),1);
    window(index3y)=min(window(index3y),0.5*(1+cos(2*pi/r*(yy(index3y)-1+r/2))));
    
%     window(index1z)=min(window(index1z),0.5*(1+cos(2*pi/r*(zz(index1z)-r/2))));
%     window(index2z)=min(window(index2z),1);
%     window(index3z)=min(window(index3z),0.5*(1+cos(2*pi/r*(zz(index3z)-1+r/2))));
    
end