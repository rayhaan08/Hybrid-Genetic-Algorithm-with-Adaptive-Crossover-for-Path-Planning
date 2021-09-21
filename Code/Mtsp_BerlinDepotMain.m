function [y y1] = Mtsp_BerlinDepotMain(xx,cities,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Berlin 52
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, ..., x6]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright - Rayhaan Iqbal (2020)
% ADAMS Lab, UB

n=51;
dep=[565.0 575.0];%THE FIRST CITY
[r c]=size(xx);
for k=1:r
    t=1;
    for i=1:m
        do=norm(dep(1,:)-cities(xx(k,t),:));
        Dist=0;
        if xx(k,n+i)==1
            t=t+1;
            Dist=0;
            de=do;
        else
            for j=1:xx(k,n+i)-1
                Dist=Dist+norm(cities(xx(k,t),:) - cities(xx(k,t+1),:));
                t=t+1;
            end
            de=norm(cities(xx(k,t),:)-dep(1,:));  %%%%%return to starting city
            t=t+1;
        end
        Tot_dist(i)=do+Dist+de;
    end
    for  i=1:m
        y1(k,i)=Tot_dist(i);
        
    end
    y(k)=max(Tot_dist);
end
y=y*-1;
end
