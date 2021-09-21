function [Opt,Opt_design]= GA_Berlin(m,runs)

% ------------------------------------------------------------------------------------------------------------------
% INPUTS
% m: No of salesman
% r: No of runs to be performed
% ------------------------------------------------------------------------------------------------------------------
% OUTPUTS
% Opt: Optimum longest tour
% Opt_design: Optimum route in two-part chromosome form
% ------------------------------------------------------------------------------------------------------------------
%
% Copyright - Rayhaan Iqbal (2020)
% ADAMS Lab, UB

fitness=@Mtsp_BerlinDepotMain;

city=[25.0 185.0; 345.0 750.0; 945.0 685.0; 845.0 655.0; 880.0 660.0; 25.0 230.0; 525.0 1000.0; 580.0 1175.0; 650.0 1130.0; 1605.0 620.0 ; 1220.0 580.0; 1465.0 200.0; 1530.0 5.0; 845.0 680.0; 725.0 370.0; 145.0 665.0; 415.0 635.0; 510.0 875.0 ;  560.0 365.0; 300.0 465.0; 520.0 585.0; 480.0 415.0; 835.0 625.0; 975.0 580.0; 1215.0 245.0; 1320.0 315.0; 1250.0 400.0; 660.0 180.0; 410.0 250.0; 420.0 555.0; 575.0 665.0; 1150.0 1160.0; 700.0 580.0; 685.0 595.0; 685.0 610.0; 770.0 610.0; 795.0 645.0; 720.0 635.0; 760.0 650.0; 475.0 960.0; 95.0 260.0; 875.0 920.0; 700.0 500.0; 555.0 815.0; 830.0 485.0; 1170.0 65.0; 830.0 610.0; 605.0 625.0; 595.0 360.0; 1340.0 725.0; 1740.0 245.0];
n=51;Npop=1000;


run_max=zeros(runs,1);run_design=zeros(runs,n+m);
for run=1:runs
    tic
    
    %%% CREATING INITIAL CITIES
    a= 1:n;
    for i=1:Npop
        b(i,randperm(n)) = a(1,:);
    end
    
    %%% INITIALLY ASSIGNING CITIES TO SALESMAN
    aa=[];
    for j=1:Npop
        remsum=n-(m-1); %%%%To ensure atleast 1 for each salesman
        for i=1:m-1
            aa(j,i)=randi([1 remsum])  ;  %%%%%1 here as minimum of 1 cities
            remsum=n-sum(aa(j,:))-(m-i-1);    %%%% (m-i-1) to account for atleast 1 city for remaining salesmen
        end
        aa(j,i+1)=n-sum(aa(j,:));
    end
    
    x=[b aa];
    [fit fit_salesman]=fitness(x,city,m);  %%%%% Here negative of objective function is called to convert to maximization problem
    xnew=x;
    
    itr_max=[];itr_avg=[];
    itrmax=1000;tolmin=0.05;tolmax=0.2;
    for itr=1:itrmax
        itr;
        count1=0;count2=0;
        sp_tol=tolmin + itr*((tolmax-tolmin)/itrmax);
        %% Roulette Wheele Selction        needs fitness to be in maximamization form
        if min(fit)==max(fit)
            break
        end
        fit_rw=((fit-min(fit))/(max(fit)-min(fit))).*(1-0)  + 0;
        prop_i=fit_rw/sum(fit_rw);
        Ni=floor(prop_i.*Npop);
        wheel=cumsum(prop_i.*360);
        while 1
            if sum(Ni)==Npop
                break
            end
            ind=find(wheel>=randi([0,360]),1);
            Ni(ind)=Ni(ind)+1;
        end
        checkingg=itr;
        %Creating the mating pool
        Mp=[];
        for i=1:Npop
            Mp=vertcat(Mp,repmat(xnew(i,:),Ni(i),1));
        end
        Mp_copy=Mp;
        
        %% CROSSOVER REAL DESIGN VARIABLE SPACE
        
        cr_prob=rand;
        if (cr_prob<0.85)
            
            count=1;seg_ind=[];seg_ind2=[];
            for i=1:Npop/2
                [r c]=size(Mp);
                cross_ind=randperm(r,2);;
                mom=Mp(cross_ind(1),:);
                dad=Mp(cross_ind(2),:);
                [Crossover Crossover_Index]=ismember(mom,Mp_copy,'rows');
                [Crossover2 Crossover_Index2]=ismember(dad,Mp_copy,'rows');
                Mp(cross_ind,:)=[];
                tol=(max(fit_salesman(Crossover_Index,:))-mean(fit_salesman(Crossover_Index,:)))/mean(fit_salesman(Crossover_Index,:));
                tol2=(max(fit_salesman(Crossover_Index2,:))-mean(fit_salesman(Crossover_Index2,:)))/mean(fit_salesman(Crossover_Index2,:));
                
                if tol<sp_tol
                    count1=count1+1;
                    %% FIRST CHILD
                    cross_pt=randi([1 n],1,2);
                    cross_pt=sort(cross_pt);
                    
                    if cross_pt(1)==1
                        UnSeg=mom(cross_pt(2)+1:n);
                    else
                        UnSeg=[mom(1:cross_pt(1)-1) mom(cross_pt(2)+1:n)];
                    end
                    if isempty(UnSeg)==0
                        
                        t=1;UnSeg_ind=[];UnSeg_copy=[];
                        for j=1:n
                            UnSeg_ind=find(UnSeg==dad(j));
                            if(UnSeg_ind>0)
                                UnSeg_copy(t)=UnSeg(UnSeg_ind);
                                t=t+1;
                            end
                        end
                        
                        size1=cross_pt(1)-1;
                        size2=n-cross_pt(2);
                        ch(count,1:size1)=UnSeg_copy(1:size1);
                        ch(count,cross_pt(1):cross_pt(2))=mom(cross_pt(1):cross_pt(2));
                        ch(count,cross_pt(2)+1:n)=UnSeg_copy(size1+1:end);
                    else
                        ch(count,1:n)=mom(1:n);
                    end
                    sales_crosspt=randi(m);
                    sales_cross=mom(n+1:n+m);
                    ch(count,(n+1:n+m))=[flip(sales_cross(1:sales_crosspt-1)) flip(sales_cross(sales_crosspt:m))];
                else
                    %% Mom cross over part 2
                    count2=count2+1;
                    t=0;seg=[];
                    for j=1:m
                        seg_ind=randi([1 mom(n+j)],1,2);    %%%possible to not select any segment from a salesmans cities
                        seg_ind=sort(seg_ind);
                        %%%%%if seg_ind==0   provision for zero segments of a salesman
                        seg_temp=mom(t+seg_ind(1):t+seg_ind(2));
                        salesman_ind(j)=length(seg_temp);
                        seg=[seg seg_temp];
                        t=t+mom(n+j);
                        
                        
                    end
                    unsaved=setdiff(mom(1:n),seg);
                    
                    if isempty(unsaved)==0
                        t=1;unsaved_ind=[];unsaved_copy=[];
                        for j=1:n
                            unsaved_ind=find(unsaved==dad(j));
                            if(unsaved_ind>0)
                                unsaved_copy(t)=unsaved(unsaved_ind);
                                t=t+1;
                            end
                        end
                        
                        ch1=[];t=1;
                        salesman_indCopy=salesman_ind;
                        for j=1:m
                            if isempty(unsaved_copy)
                                ch1=[ch1 seg(t:end)];
                                break
                            end
                            
                            ch1=[ch1 seg(t:(t+salesman_indCopy(j)-1))];
                            uscr(j)=randi([1 length(unsaved_copy)]);
                            ch1=[ch1 unsaved_copy(1:uscr(j))];
                            unsaved_copy(1:uscr(j))=[];
                            t=t+salesman_indCopy(j);
                            salesman_ind(j)=salesman_indCopy(j)+uscr(j);
                            
                        end
                        ch1=[ch1 unsaved_copy];
                        salesman_ind(m)=salesman_ind(m)+length(unsaved_copy);
                    else
                        ch1=seg;
                    end
                    ch(count,:)=[ch1 salesman_ind];
                end
                
                if tol2<sp_tol
                    %% SECOND CHILD
                    count1=count1+1;
                    cross_pt2=randi([1 n],1,2);
                    cross_pt2=sort(cross_pt2);
                    if cross_pt2(1)==1
                        UnSeg2=dad(cross_pt2(2)+1:n);
                    else
                        UnSeg2=[dad(1:cross_pt2(1)-1) dad(cross_pt2(2)+1:n)];
                    end
                    if isempty(UnSeg2)==0
                        
                        t=1;UnSeg_ind2=[];UnSeg_copy2=[];
                        for j=1:n
                            UnSeg_ind2=find(UnSeg2==mom(j));
                            if(UnSeg_ind2>0)
                                UnSeg_copy2(t)=UnSeg2(UnSeg_ind2);
                                t=t+1;
                            end
                        end
                        
                        size1=cross_pt2(1)-1;
                        size2=n-cross_pt2(2);
                        ch(count+1,1:size1)=UnSeg_copy2(1:size1);
                        ch(count+1,cross_pt2(1):cross_pt2(2))=dad(cross_pt2(1):cross_pt2(2));
                        ch(count+1,cross_pt2(2)+1:n)=UnSeg_copy2(size1+1:end);
                    else
                        ch(count+1,1:n)=dad(1:n);
                    end
                    sales_crosspt2=randi(m);
                    sales_cross2=dad(n+1:n+m);
                    ch(count+1,(n+1:n+m))=[flip(sales_cross2(1:sales_crosspt2-1)) flip(sales_cross2(sales_crosspt2:m))];
                else
                    %% dad crossover part 2
                    count2=count2+1;
                    seg2=[];t2=0;
                    for j=1:m
                        seg_ind2=randi([1 dad(n+j)],1,2);    %%%possible to not select any segment from a salesmans cities
                        seg_ind2=sort(seg_ind2);
                        %%%%%if seg_ind==0   provision for zero segments of a salesman
                        seg_temp2=dad(t2+seg_ind2(1):t2+seg_ind2(2));
                        salesman_ind2(j)=length(seg_temp2);
                        seg2=[seg2 seg_temp2];
                        t2=t2+dad(n+j);
                    end
                    unsaved2=setdiff(dad(1:n),seg2);
                    
                    if isempty(unsaved2)==0      %%% if not empty
                        t2=1;unsaved_ind2=[];unsaved_copy2=[];
                        for j=1:n
                            unsaved_ind2=find(unsaved2==mom(j));
                            if(unsaved_ind2>0)
                                unsaved_copy2(t2)=unsaved2(unsaved_ind2);
                                t2=t2+1;
                            end
                        end
                        unsc2=unsaved_copy2;
                        
                        ch2=[];t2=1;
                        salesman_indCopy2=salesman_ind2;
                        for j=1:m
                            if isempty(unsaved_copy2)
                                ch2=[ch2 seg2(t2:end)];
                                break
                            end
                            ch2=[ch2 seg2(t2:(t2+salesman_indCopy2(j)-1))];
                            uscr2(j)=randi([1 length(unsaved_copy2)]);
                            ch2=[ch2 unsaved_copy2(1:uscr2(j))];
                            unsaved_copy2(1:uscr2(j))=[];
                            t2=t2+salesman_indCopy2(j);
                            salesman_ind2(j)=salesman_indCopy2(j)+uscr2(j);
                        end
                        ch2=[ch2 unsaved_copy2];
                        salesman_ind2(m)=salesman_ind2(m)+length(unsaved_copy2);
                    else
                        ch2=seg2;
                    end
                    ch(count+1,:)=[ch2 salesman_ind2];
                end
                count=count+2;
            end
            %% if no crossover
        else
            ch=Mp_copy;
            xnew=ch;
        end
        
        
        %% Mutation
        mut_prob=rand;
        if mut_prob<0.9
            for k=1:Npop
                mut=randi([1 n],1,2);
                temp=ch(k,mut(1));
                ch(k,mut(1))=ch(k,mut(2));
                ch(k,mut(2))=temp;
                
                mut=randi([1 m],1,2);
                temp=ch(k,n+mut(1));
                ch(k,n+mut(1))=ch(k,n+mut(2));
                ch(k,n+mut(2))=temp;
            end
        end
        
        
        
        
        %% ELITISM
        if itr>20
            [fit_ch fit_SalesCh]=fitness(ch,city,m);
            fit_temp=[fit_ch fit];
            x_temp=[ch; xnew];
            fit_SalesTemp=[fit_SalesCh; fit_salesman];
            [fit_temp indxs]=sort(fit_temp,'descend');
            
            fit_salesman=fit_SalesTemp(indxs(1:Npop),:);
            xnew=x_temp(indxs(1:Npop),:);
            fit=fit_temp(1,1:Npop);
            %xnew=ch;%%%%%%%% IMPLEMENT NSGA
        else
            [fit fit_salesman]=fitness(ch,city,m);
            xnew=ch;
        end
        
        itr_count1(itr)=count1;
        itr_count2(itr)=count2;
        
        %xnew=ch;
        itr_avg(itr)=-mean(fit);
        [itr_max(itr) dsign]=min(-fit); %%%%As the obj function is called as negative so max prob
        itr_design(itr,:)=xnew(dsign,:);
        itr_totdis(itr)=sum(fit_salesman(dsign,:));
    end
    [run_max(run) dsign]=min(itr_max);
    run_design(run,:)=itr_design(dsign,:);
    run_totdis(run)=itr_totdis(dsign);
    
    
    
    time_run(run)=toc;
    time_runItr(run)=itr;
end

[disp disp_ind]=min(run_max);
plotting=run_design(disp_ind,52:51+m);
a=1;
figure
temp=[];xplot=[];yplot=[];
for i=1:m
    temp=[];xplot=[];yplot=[];
    temp=run_design(disp_ind,a:a+plotting(i)-1);
    a=a+plotting(i);
    xplot=[565.0 city(temp,1)' 565.0];
    yplot=[575.0 city(temp,2)' 575.0];
    plot(xplot,yplot)
    hold on
    
end

plot(city(:,1),city(:,2),'o')
hold on
plot(565,576,'x')
xlabel('X-Coordinate')
ylabel('Y-Coordinate')
title('Berlin52')


%CONVERGENCE HISTORY
 figure
plot(1:length(itr_max),itr_max)
ylabel('Iteration Best')
xlabel('Iteration')
title({'Berlin52','Convergence History','Iteration Best vs Iteration No.'})

figure
plot(1:length(itr_avg),itr_avg)
ylabel('Iteration Average')
xlabel('Iteration')
title({'Berlin52','Convergence History','Iteration Average vs Iteration No.'})

Opt_design=run_design(disp_ind,:);
Opt=min(run_max);
 end



