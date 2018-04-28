function [centroid,b]=segmentClathrin(im)
    D_new=zeros(size(im));
    im=imadjust(im);
    img=imgaussfilt(im,1);
    %im2=medfilt2(img);
    im2=im2double(img);
    h=fspecial('log',12,1.5);
    imf=imfilter(im2,h,'same');
    imf1=imcomplement(imf);
    l=multithresh(imf1);
    bw=imbinarize(imf1,l);
    se=strel('disk',1);
%     bw1=imclose(bw,se);
    newbw=imdilate(bw,se);
%     [gmag,gdir]=imgradient(imf1);
%     [D,L]=bwdist(~bw1,'euclidean');
%     D=-D;
%     D(~bw1)= Inf;
%     gmax=max(gmag(:));
%     gmin=min(gmag(:));
%     for i=1:size(gmag,1)
%         for j=1:size(gmag,2)
%             D_new(i,j)=D(i,j)*exp(1-(gmag(i,j)-gmin)/(gmax-gmin));
%         end
%     end
%     
%     DL=watershed(D_new);
%     DL(~bw1)=0;
%     t=adaptthresh(DL,1);
%     newbw=imbinarize(DL,t);
%     thresh1=0.3;thresh2=25;
%     thresh3=0.5;thresh4=33;
    [b,lab,n,a]=bwboundaries(newbw);
    cc=bwconncomp(newbw);
    stats=regionprops(cc,'Area','Centroid');
    centroid=cat(1,stats.Centroid);
end