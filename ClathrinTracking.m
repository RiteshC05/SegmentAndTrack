clear all;
imt='loci1.tif';
info=imfinfo(imt);
num=numel(info);
T=1;
accel=0;
R=[0.1 0;0 0.1];
Q=eye(4);
P=Q;
A=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];
B=[0.5*T^2 ;0.5*T^2;T;T];
C=[1 0 0 0;0 1 0 0];
im=imread('loci1.tif');
[centroids,b]=segmentClathrin(im);
 X=centroids(:,1);
 Y=centroids(:,2);
 bbox=cell(600,1500);
 bbox(1,1:size(b,1))=b;
 initial_est=[X(:)';Y(:)';(zeros(size(X)))';(zeros(size(Y)))']; 
 num_detections=size(X,1);
 num_tracks=num_detections;
 x_coords=nan(600,1500);
 y_coords=nan(600,1500);
 tracks_col=2;
 tracks=[X Y];
 penalty_track=zeros(3000,1);
for k=1:num
    im = imread('loci1.tif', k, 'Info', info);
    [centroids,b]=segmentClathrin(im);
    X=centroids(:,1);
    Y=centroids(:,2);
    present_obsv=[X Y];%getting new observations
    intial_est=A*initial_est + B*accel;%finding new estimate
    P=A*P*A' +Q;%covariance matrix
    K=P*C'*inv(C*P*C'+R);%kalman gain
    dist_map=pdist([initial_est(1:2,:)';present_obsv]);%calculating distance map
    distmap=squareform(dist_map);
    distmap=distmap(1:num_tracks,num_tracks+1:end);
    [assignment,cost]=assignmentoptimal(distmap);%hungarian algorithm
    assgn_length=length(assignment);
    detect_length=size(X,1);
    newtemp=[]; 
    for i=1:assgn_length
       if (assignment(i)>0 && distmap(i,assignment(i))>50)
           assignment(i)=0;
       end
    end
    newbbox=bbox(k,1:num_tracks);
    %newtemp=present_obsv(assignment,1:2);
    %for i=1:assgn_length
     %   temp=cat(1,temp,newtemp);
    %end
    %initial_est=new_est;
    j=1;
    for i=1:assgn_length
        if(assignment(i)>0)
        initial_est(:,j)=initial_est(:,j)+K*(present_obsv(assignment(i),:)'-C*initial_est(:,j));
        newbbox(1,j)=b(assignment(i),1);
        end
        j=j+1;
    end

    P=(eye(4)-K*C)*P;
    if (ismember(0,assignment))
        i=1:size(present_obsv,1);
         ind=~ismember(i,assignment);
         new_track=[];
         new_track=present_obsv(ind,:)';
         new_box=b(ind,1);
         new_track=cat(1,new_track,zeros(2,size(new_track,2)));
         initial_est=cat(2,initial_est,new_track);
         num_tracks=size(initial_est,2); 
         newbbox=cat(2,newbbox,new_box');
    end
    x_coords(k,1:num_tracks)=initial_est(1,1:num_tracks);
    y_coords(k,1:num_tracks)=initial_est(2,1:num_tracks);
    bbox(k,1:num_tracks)=newbbox(1,1:num_tracks);
        
    %if (detect_length>assgn_length)
        % i=1:size(present_obsv,1);
        % ind=~ismember(i,assignment);
        % new_track=[];
        % new_track=present_obsv(ind,:)';
       %  new_track=cat(1,new_track,zeros(2,size(new_track,2)));
      %   initial_est=cat(2,initial_est,new_track);
     %    num_tracks=size(initial_est,2);  
    %end
     
     for i=1:assgn_length
         if assignment(i)==0
             penalty_track(i)=penalty_track(i)+1;
         end
     end
     
     penalty_ind=find(penalty_track>8);
     for i=1:length(penalty_ind)
         initial_est(:,penalty_ind(i))=NaN;
     end
             
             
     
    %norm_track=present_obsv(ind,:);
    %new_track1=cat(1,new_track,zeros(size(new_track)));
    %diff_col=size(tracks,2)-size(norm_track,2);
    %new_track2=cat(2,zeros(size(norm_track,1),diff_col),norm_track);
    %initial_est=cat(2,updated_est,new_track1);
    %tracks=cat(1,tracks,new_track2);
    %tracks_col=size(tracks,2);
  
end
for i=1:length(penalty_track)
    if (penalty_track(i)>250)
        x_coords(i,:)=nan;
        y_coords(i,:)=nan;
    end
end
% dist1=[];
% for j=1:num_tracks
%     for i=1:600
%         if (~isnan(x_coords(i,j)))
%             dist1(i,j)=sqrt((1-x_coords(i,j))^2 + (1-y_coords(i,j))^2);
%         else
%             dist1(i,j)=0;
%         end
%     end
% end
% dist1=sum(dist1);
% dist1=dist1./1000;
% med=median(dist1);
% j=1;
% v=VideoWriter('Trackingfinal14.avi');
% v.FrameRate=12;
% open(v);
cmap=hsv(255);
for k=1:600
     %colors={'b','r','g','y','k'};
     im = imread('loci1.tif', k, 'Info', info);
     [centroids]=segmentClathrin(im);
     X=centroids(:,1);
     Y=centroids(:,2);
     figure(1);
     im=imadjust(im);
     imagesc(im);
     hold on
     for i=1:num_tracks
%           idx=1+round(rand()*254);
         if (~isnan(x_coords(k,i)))
             bound=bbox{k,i};
             if (~isempty(bound)) 
                 plot(bound(:,2),bound(:,1),'c','LineWidth',2);
             
%              if (dist(i)>med)
%                  plot(x_coords(k,i),y_coords(k,i),'g*');
%              else
%                  plot(x_coords(k,i),y_coords(k,i),'c*');
             end
             text(x_coords(k,i)+0.2, y_coords(k,i)+0.2,num2str(i),'Color','white','FontSize',10);
                 
          end
     end
     hold off
%      writeVideo(v,getframe);
   %  j=j+1;
end

% close(v)