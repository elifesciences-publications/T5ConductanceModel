function optimize_model(task_id)

% example function call:
% optimize_model('42')
% will use 42 as a seed for the random number generator and 
% create a folder named '42' with optimization results for all
% recorded T5 cells. The results can be loaded and visualized with the code
% available at
% https://doi.org/10.25378/janelia.c.4771805.v1

attempts_num = 10;
err_cutoff = 3.5;
basedir ='./';

cd(basedir)
newdir = strcat(basedir,task_id);
mkdir(newdir)
rng(str2double(task_id));

lb_gauss = [ .1,.1,.1,.1,.1,   0,0,-5,-5,0,0];
ub_gauss = [ 40,40,40,40,40,   10,10,5,10,10,10];

for cell_id =  1:17
    load([basedir,'data_cell_',num2str(cell_id),'_spfr.mat'],'d','p')
    for n_attemps = 1:attempts_num
        par0 = lb_gauss + rand(size(ub_gauss)).*(ub_gauss-lb_gauss);
        q = optimize_t5(d,p,par0,lb_gauss,ub_gauss);
        vm_spfr = vm_bars_and_gratings(q.param,0,p);
        err_spfr = mean( (vm_spfr - d.vm_all_sb).^2 );
        if(q.resnorm <  err_cutoff)
            break
        end
    end
    load([basedir,'data_cell_',num2str(cell_id),'_all.mat'],'d','p') 
    vm_all = vm_bars_and_gratings(q.param,0,p);
    if(ismember(cell_id, [2 5 6 7 9 10 13 14 15 16 17]))
        err_all = mean( (vm_all -...
            [d.vm_all_sb; d.vm_all_mb; d.vm_all_mm;...
            d.vm_all_mg; d.vm_all_sg]).^2 );
    else
        err_all = mean( (vm_all -...
            [d.vm_all_sb; d.vm_all_mb; d.vm_all_mm]).^2 );
    end
    save([newdir,'/result_cell_',num2str(cell_id),'.mat'],...
        'par0','q','vm_all','err_spfr','err_all')
end

function q = optimize_t5(d,p,par0,lb,ub)
vm_all = d.vm_all_sb; 
ind = false(size(vm_all));
for i = 1:length(d.ind_ref1_sb)
    wid = p.width_sb(i);
    if(wid==2)
        ind(d.ind_ref1_sb(i):d.ind_ref2_sb(i)) = true;
    end
end

fun = @(ppar) mean( (extract_rise_decay_vm(ppar, p, ind) -...
    vm_all(ind) ).^2 );  
options = optimoptions('fmincon','Display','iter');
[param,resnorm,residual,exitflag,output] =...
    fmincon(fun,par0,[],[],[],[],lb,ub,[],options);
   
q.param = param;
q.resnorm = resnorm;
q.residual = residual;
q.exitflag = exitflag;
q.output = output;

function vm_ind = extract_rise_decay_vm(par, p, ind)
tmp = vm_bars_and_gratings(par, 0, p);
vm_ind = tmp(ind);

function y = vm_bars_and_gratings(param,~,q)

p.xx = q.xx;
p.tau = param(1)*10;
p.tde = param(2)*10;
p.tdi = param(3)*10;
p.tre = param(4)*10;
p.tri = param(5)*10;

ae = param(6);
ai = param(7);
mue = param(8);
mui = param(9);
sige = param(10);
sigi = param(11);
p.xae = ae * exp(-(p.xx-mue).^2 / (2 *sige^2));
p.xai = ai * exp(-(p.xx-mui).^2 / (2 *sigi^2));

p.vi =  -74;
p.vl = -65; 
p.ve =  0;
p.toff = 25;
p.stimmin = p.toff;

y = [];
k0 = 1; k1 = 1; k2 = 1; k3 = 1; k4 = 1;
for i = 1:length(q.protocol)
    incond = [0,0,0,0];
    
    if(q.protocol(i)==0) % single bars
        p.tstimmin = p.toff;
        p.tstimmax = q.duration_sb(k0) + p.toff;
        p.location_sb = q.location_sb(k0);       
        p.width_sb = q.width_sb(k0);
        [~,y0] = ode45(@run_vm_sb, q.t{i}, incond, [],p);
        vm = (p.vl + p.ve*y0(:,2) + p.vi*y0(:,4))./(1 + y0(:,2) + y0(:,4));
        y = [y; vm - p.vl];
        k0 = k0 + 1;
    end
    
    if(q.protocol(i)==3) % moving bars
        p.pos = q.pos_mb;
        p.flag_dir = q.direction_mb(k3);
        p.currwidth = q.width_mb(k3);
        p.currdur = q.duration_mb(k3);        
        [~,y0] = ode45(@run_vm_mb, q.t{i}, incond, [],p);
        vm = (p.vl + p.ve*y0(:,2) + p.vi*y0(:,4))./(1 + y0(:,2) + y0(:,4));
        y = [y; vm - p.vl];
        k3 = k3 + 1;
    end
    
    if(q.protocol(i)==4) % minimal motion
        p.fb_pos = q.fb_pos_mm(k4);
        p.sb_pos = q.sb_pos_mm(k4);
        p.width = q.width_mm(k4);
        p.fbton = q.fbton_mm + p.toff;
        p.fbtoff = q.fbtoff_mm(k4) + p.toff;
        p.sbton = q.sbton_mm(k4) + p.toff;
        p.sbtoff = q.sbtoff_mm(k4) + p.toff;
        [~,y0] = ode45(@run_vm_mm, q.t{i}, incond, [],p);
        vm = (p.vl + p.ve*y0(:,2) + p.vi*y0(:,4))./(1 + y0(:,2) + y0(:,4));
        y = [y; vm - p.vl];
        k4 = k4 + 1;
    end
    
    if(q.protocol(i)==1) % moving grating
        p.PosD_mg = q.PosD_mg;
        p.currdur = q.stimdur_mg(k1);
        p.phaseMat = q.phase_mg(k1,:);
        [~,y0] = ode45(@run_vm_mg, q.t{i}, incond, [],p);  
        vm = (p.vl + p.ve*y0(:,2) + p.vi*y0(:,4))./(1 + y0(:,2) + y0(:,4));
        y = [y; vm - p.vl];
        k1 = k1 + 1;
    end
    
    if(q.protocol(i)==2) % static grating
        p.tstimmin = p.toff;
        p.tstimmax = p.toff + q.stimdur_sg(k2);
        p.posD_sg1 = q.posD_sg{k2};
        [~,y0] = ode45(@run_vm_sg, q.t{i}, incond, [],p);
        vm = (p.vl + p.ve*y0(:,2) + p.vi*y0(:,4))./(1 + y0(:,2) + y0(:,4));
        y = [y; vm - p.vl];
        k2 = k2 + 1;
    end
    
    
end

function dydt = run_vm_sb(t,y,p) % single bars

stim = 0;
if(t>=p.tstimmin && t<=p.tstimmax)
    stim = 1;
end
ae = 0; ai = 0;
j = find(p.xx==p.location_sb);        
for k = 0:(p.width_sb-1)
    ae = ae + p.xae(j-k);
    ai = ai + p.xai(j-k);
end
ae = ae*stim;
ai = ai*stim;
dydt = rhs_vm(ae,ai,y,p);

function dydt = run_vm_mm(t,y,p) % minimal motion

stim1 = 0;
if(t>=p.fbton && t<=p.fbtoff)
    stim1 = 1;
end
stim2 = 0;
if(t>=p.sbton && t<=p.sbtoff)
    stim2 = 1;
end

xflag = zeros(size(p.xae));
ae = 0; ai = 0;
if(stim1==1)
    j = find(p.xx==p.fb_pos); 
    for k = 0:(p.width-1)
        ae = ae + p.xae(j-k);
        ai = ai + p.xai(j-k);
        xflag(j-k) = 1;
    end
end

if(stim2==1)
    j = find(p.xx==p.sb_pos); 
    for k = 0:(p.width-1)
        if(xflag(j-k)~=1)
            ae = ae + p.xae(j-k);
            ai = ai + p.xai(j-k);
        end
    end
end

dydt = rhs_vm(ae,ai,y,p);


function dydt = run_vm_mb(t,y,p) %moving bar

istim = floor((t-p.toff)/(  p.currdur));
ae = 0; ai = 0;  

eff_pos = p.pos;
if(p.flag_dir==1)
    for i=1:(p.currwidth-1)
        eff_pos = [eff_pos p.pos(end)+i] ;
    end
end
if(p.flag_dir==0)
    for i=1:(p.currwidth-1)
        eff_pos = [p.pos(1)-i eff_pos] ;
    end
    eff_pos = fliplr(eff_pos);
end 

if(t>p.toff && istim<length(eff_pos))
    lim_width = (p.currwidth-1);
    if(p.flag_dir==1)
        sloc = find(p.xx==eff_pos(istim+1));
        if(istim==0); lim_width = 0; end
        if(istim==1); lim_width = min(1,p.currwidth-1); end
        if(istim==2); lim_width = min(2,p.currwidth-1); end
        if(istim>=length(p.pos))
            sloc = find(p.xx==p.pos(end));
            if(istim==length(p.pos)); lim_width = p.currwidth-2; end
            if(istim==length(p.pos)+1); lim_width = p.currwidth-3; end
            if(istim==length(p.pos)+2); lim_width = p.currwidth-4; end
        end
        sign_fact=-1;       
    end
    if(p.flag_dir==0)
        sloc = find(p.xx==eff_pos(istim+1));
        if(istim==0); lim_width = 0; end
        if(istim==1); lim_width = min(1,p.currwidth-1); end
        if(istim==2); lim_width = min(2,p.currwidth-1); end
        if(istim>=length(p.pos))
            sloc = find(p.xx==p.pos(1));
            if(istim==length(p.pos)); lim_width = p.currwidth-2; end
            if(istim==length(p.pos)+1); lim_width = p.currwidth-3; end
            if(istim==length(p.pos)+2); lim_width = p.currwidth-4; end
        end
        sign_fact=1;    
    end
    for k = 0:lim_width
        ae = ae + p.xae(sloc+ sign_fact*k);
        ai = ai + p.xai(sloc+ sign_fact*k);
    end
end
dydt = rhs_vm(ae,ai,y,p);


function dydt = run_vm_mg(t,y,p) % moving grating

istim = floor((t-p.toff)/(  p.currdur));
ae = 0; ai = 0;  
if(t>p.toff && istim<=3*8)
    pos_ind = p.phaseMat( mod(istim,8) + 1 );
    pos_d = p.PosD_mg{pos_ind};
    for ii = 1:length(pos_d)
        ae = ae + p.xae( p.xx==pos_d(ii) );
        ai = ai + p.xai( p.xx==pos_d(ii) );
    end
   
end
dydt = rhs_vm(ae,ai,y,p);


function dydt = run_vm_sg(t,y,p) %static grating

stim = 0;
ae = 0; ai = 0;  
for ii = 1:length(p.posD_sg1)
    ae = ae + p.xae( p.xx==p.posD_sg1(ii) );
    ai = ai + p.xai( p.xx==p.posD_sg1(ii) );
end
if(t>=p.tstimmin && t<=p.tstimmax)
    stim = 1;
end
ae = ae * stim;
ai = ai * stim;
dydt = rhs_vm(ae,ai,y,p);


function dydt = rhs_vm(ae,ai,y,p)

he = y(1);
ge = y(2);
hi = y(3);
gi = y(4);

dhedt = ( -he + ae ) / p.tre;
dgedt = ( -ge + he ) / p.tde;
dhidt = ( -hi + ai ) / p.tri;
dgidt = ( -gi + hi ) / p.tdi;
dydt = [dhedt; dgedt; dhidt; dgidt];


