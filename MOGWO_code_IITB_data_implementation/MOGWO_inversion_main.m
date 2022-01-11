% Multi-Objective grey wolf optimiser part adapted from Mirjalili et al.(2016), Mirjalili(2021)
% Rayleigh wave phase velocity dispersion curve generation inspired from Herrmann(Computer Program in Seismology, 2013) and Lehujeur et al.(2018)
% Surface wave processing inspired from Vantassel and Cox(2021a,b)
% If you use the following Rayleigh wave phase velocity joint inversion using multi-objective grey wolf optimiser codes in your research please cite Vashisth et al.(2021) 
% The following code was executed in GNU Octave and this sample does not require GPU to run

clear

[status,cmdout]=unix ("python observed_data_read1.py");
vel_phase_fundamental=load('observed_fund_phase_velocity.txt'); % reading fundamental mode phase velocities from file

[status,cmdout]=unix ("python observed_data_read2.py");
vel_phase_first=load('observed_first_phase_velocity.txt'); % reading first mode phase velocities from file


rand('state',1); %Seeding random number generator

nn=5; %number of layers in the subsurface model


zlb=[1, 4, 7, 14]; % lower bounds for depths corresponding to layer top(m)
zub=[2, 6, 9, 16]; % upper bounds for depths corresponding to layer top(m)


vplb=[500,    1500,   2500,   3000,   3500]; % lower bounds for P-wave velocities(m/s)
vpub=[1200,   2500,   3500,   4000,   4500]; % upper bounds for P-wave velocities(m/s)


vslb=[85,   200,   500,   600,   900]; % lower bounds for S-wave velocities(m/s)
vsub=[200,  400,   700,   800,   1250]; % upper bounds for S-wave velocities(m/s)


rolb=[1.01,  1.6,   2.0,   2.3,   2.4];  %if you want to invert for density(g/cm^3) along with Vs, Vp, and z. If not, comment this line.   
roub=[1.15,  1.9,   2.3,   2.6,   2.7];  %if you want to invert for density(g/cm^3) along with Vs, Vp, and z. If not, comment this line.

lb=[zlb,vplb,vslb,rolb]; %if not inverting for density, lb=[zlb,vplb,vslb]
ub=[zub,vpub,vsub,roub]; %if not inverting for density, ub=[zub,vpub,vsub]

nVar=(4*nn)-1;  %total number of unknowns. (3*nn)-1 if not inverting for density

Positions_initial=rand(1,nVar).*(ub-lb)+lb;

fobj = @(x) joint_inv_obj_fun(x,vel_phase_fundamental,vel_phase_first);

VarSize=[1 nVar];

GreyWolves_num=50; % Total no. of search agents (tuning parameter)
MaxIt=50;  % Maximum Number of Iterations (tuning parameter)
Archive_size=50;   % Repository Size

alpha=0.1;  % Grid Inflation Parameter
nGrid=10;   % Number of Grids per each Dimension
beta=4;     % Leader Selection Pressure Parameter
gamma=2;    


empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Dominated=false;
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];
    
GreyWolves=repmat(empty_particle,GreyWolves_num,1);


for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    GreyWolves(i).Position=zeros(1,nVar);
    for j=1:nVar
        GreyWolves(i).Position(1,j)=unifrnd(lb(j),ub(j),1);
    end
    GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end


npop=numel(GreyWolves);
for i=1:npop
    GreyWolves(i).Dominated=false;
    for j=1:i-1
        if ~GreyWolves(j).Dominated
            if all(GreyWolves(i).Cost<=GreyWolves(j).Cost) && any(GreyWolves(i).Cost<GreyWolves(j).Cost)
                GreyWolves(j).Dominated=true;
            elseif all(GreyWolves(j).Cost<=GreyWolves(i).Cost) && any(GreyWolves(j).Cost<GreyWolves(i).Cost)
                GreyWolves(i).Dominated=true;
                break;
            end
        end
    end
end


ND_GreyWolves=~[GreyWolves.Dominated];
Archive=GreyWolves(ND_GreyWolves); % Getting Non-dominated solutions

Archive_costs=reshape([Archive.Cost],numel(Archive(1).Cost),[]);

%% Creating Hypercubes
empty_grid.Lower=[];
empty_grid.Upper=[];
G=repmat(empty_grid,size(Archive_costs,1),1);

for j=1:size(Archive_costs,1)
    min_cj=min(Archive_costs(j,:));
    max_cj=max(Archive_costs(j,:));

    dcj=alpha*(max_cj-min_cj);

    min_cj=min_cj-dcj;
    max_cj=max_cj+dcj;

    gx=linspace(min_cj,max_cj,nGrid-1);

    G(j).Lower=[-inf gx];
    G(j).Upper=[gx inf];
end



for i=1:numel(Archive)
    %%Getting grid index
    str=['sub2ind(' mat2str(ones(1,numel(Archive(i).Cost))*numel(G(1).Upper))];

    Archive(i).GridSubIndex=zeros(1,numel(Archive(i).Cost));
    for j=1:numel(Archive(i).Cost)
        
        ij=find(Archive(i).Cost(j)<G(j).Upper,1,'first');
        
        Archive(i).GridSubIndex(j)=ij;
        
        str=[str ',' num2str(ij)];
    end
    
    str=[str ');'];
    
    Archive(i).GridIndex=eval(str);
end

% MOGWO main loop

for ii=1:MaxIt
    a=2-ii*((2)/MaxIt);  
    for i=1:GreyWolves_num
        
        clear hyp2
        clear hyp3
        
        Delta=Leader_Selection(Archive,beta);
        Beta=Leader_Selection(Archive,beta);
        Alpha=Leader_Selection(Archive,beta);
        
        % Ifnumber of solutions is less than three in the least crowded hypercube then second least crowded hypercube is selected
        if size(Archive,1)>1
            cc=0;
            for newid=1:size(Archive,1)
                if sum(Delta.Position~=Archive(newid).Position)~=0
                    cc=cc+1;
                    hyp2(cc,1)=Archive(newid);
                end
            end
            Beta=Leader_Selection(hyp2,beta);
        end
        
        % If second least crowded hypercube has only one solution then the delta leader is chosen from the third least crowded hypercube
        if size(Archive,1)>2
            cc=0;
            for newid=1:size(hyp2,1)
                if sum(Beta.Position~=hyp2(newid).Position)~=0
                    cc=cc+1;
                    hyp3(cc,1)=hyp2(newid);
                end
            end
            Alpha=Leader_Selection(hyp3,beta);
        end
        
        c=2.*rand(1, nVar); 
        D=abs(c.*Delta.Position-GreyWolves(i).Position);
        A=2.*a.*rand(1, nVar)-a;
        M1=Delta.Position-A.*abs(D);
        
        c=2.*rand(1, nVar);
        D=abs(c.*Beta.Position-GreyWolves(i).Position);
        A=2.*a.*rand()-a;
        M2=Beta.Position-A.*abs(D);
        
        c=2.*rand(1, nVar);
        D=abs(c.*Alpha.Position-GreyWolves(i).Position);
        A=2.*a.*rand()-a;
        M3=Alpha.Position-A.*abs(D);
        
        % Model Update
        GreyWolves(i).Position=(M1+M2+M3)./3;
        
        % Boundary check
        GreyWolves(i).Position=min(max(GreyWolves(i).Position,lb),ub);
        
        GreyWolves(i).Cost=fobj(GreyWolves(i).Position')';
    end
    
    npop=numel(GreyWolves);
    for i=1:npop
        GreyWolves(i).Dominated=false;
        for j=1:i-1
            if ~GreyWolves(j).Dominated
                if all(GreyWolves(i).Cost<=GreyWolves(j).Cost) && any(GreyWolves(i).Cost<GreyWolves(j).Cost)
                    GreyWolves(j).Dominated=true;
                elseif all(GreyWolves(j).Cost<=GreyWolves(i).Cost) && any(GreyWolves(j).Cost<GreyWolves(i).Cost)
                    GreyWolves(i).Dominated=true;
                    break;
                end
            end
        end
    end

    ND_GreyWolves=~[GreyWolves.Dominated];
    non_dominated_wolves=GreyWolves(ND_GreyWolves); 
    
    Archive=[Archive
        non_dominated_wolves];
    
    npop=numel(Archive);
    for i=1:npop
        Archive(i).Dominated=false;
        for j=1:i-1
            if ~Archive(j).Dominated
                if all(Archive(i).Cost<=Archive(j).Cost) && any(Archive(i).Cost<Archive(j).Cost)
                    Archive(j).Dominated=true;
                elseif all(Archive(j).Cost<=Archive(i).Cost) && any(Archive(j).Cost<Archive(i).Cost)
                    Archive(i).Dominated=true;
                    break;
                end
            end
        end
    end

    ND_Archive=~[Archive.Dominated];
    Archive=Archive(ND_Archive);

    for i=1:numel(Archive)
        %%Getting grid index
        str=['sub2ind(' mat2str(ones(1,numel(Archive(i).Cost))*numel(G(1).Upper))];

        Archive(i).GridSubIndex=zeros(1,numel(Archive(i).Cost));
        for j=1:numel(Archive(i).Cost)
        
            ij=find(Archive(i).Cost(j)<G(j).Upper,1,'first');
        
            Archive(i).GridSubIndex(j)=ij;
        
            str=[str ',' num2str(ij)];
        end
    
        str=[str ');'];
    
        Archive(i).GridIndex=eval(str);
    end
    
    if numel(Archive)>Archive_size
        for k=1:numel(Archive)-Archive_size
            GridIndices=[Archive.GridIndex];
            occ_cell_index=unique(GridIndices);
            occ_cell_member_count=zeros(size(occ_cell_index));

            m=numel(occ_cell_index);
            for k=1:m
                occ_cell_member_count(k)=sum(GridIndices==occ_cell_index(k));
            end

            p=occ_cell_member_count.^gamma;
            p=p/sum(p);

            selected_cell_index=occ_cell_index(find((rand)<=cumsum(p),1,'first'));

            selected_cell_members=find(GridIndices==selected_cell_index);

            selected_member_index=randi([1 numel(selected_cell_members)]);

            Archive=[Archive(1:selected_cell_members(selected_member_index)-1); Archive(selected_cell_members(selected_member_index)+1:end)];
        end
        
        Archive_costs=reshape([Archive.Cost],numel(Archive(1).Cost),[]);
   
        empty_grid.Lower=[];
        empty_grid.Upper=[];
        G=repmat(empty_grid,size(Archive_costs,1),1);

        for j=1:size(Archive_costs,1)
            min_cj=min(Archive_costs(j,:));
            max_cj=max(Archive_costs(j,:));

            dcj=alpha*(max_cj-min_cj);

            min_cj=min_cj-dcj;
            max_cj=max_cj+dcj;

            gx=linspace(min_cj,max_cj,nGrid-1);

            G(j).Lower=[-inf gx];
            G(j).Upper=[gx inf];
        end
        
    end
    
    disp(['In iteration ' num2str(ii) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    save results
    
    % Results
    
    costs=reshape([GreyWolves.Cost],numel(GreyWolves(1).Cost),[]);

    Archive_costs=reshape([Archive.Cost],numel(Archive(1).Cost),[]);
    
    
end

rand('state','reset');

i=1;
while i<=length(Archive)
mi=getfield(Archive,{i},'Position')./1000;
mi(length(mi)-nn+1:length(mi))=mi(length(mi)-nn+1:length(mi)).*1000; %density should be in g/cm^3
name=sprintf('final_model%d.txt',i);
dlmwrite(name,mi,'\n');
i=i+1;
end


