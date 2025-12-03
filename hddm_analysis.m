%This code is used once the RLDDM is ran and you have the 4 traces

%temp_tbl = full_tbl;
%open setshift_tbl.mat
protocol_tbl = full_tbl(full_tbl.Protocol=='cog_SET SHIFT_Training',:); %change based on what sessions/schedules you're running
%temp_tbl = protocol_tbl(~isnan(protocol_tbl.frontChoice),:); %this removes omissions, probably should include them.
temp_tbl = protocol_tbl;
hddm_tbl = table();
[C,~,ib] = unique(temp_tbl.Rat);
hddm_tbl.rt = temp_tbl.RT;
hddm_tbl.subj_idx = ib;
hddm_tbl.split_by = str2double(temp_tbl.Session);
hddm_tbl.sex = strcmp(temp_tbl.Sex, 'F');
hddm_tbl.light = 2 - temp_tbl.light;
hddm_tbl.feedback = temp_tbl.performance;
hddm_tbl.response = 2-temp_tbl.frontChoice;
hddm_tbl.response(isnan(temp_tbl.frontChoice)) = -1;
hddm_tbl.rt(isnan(temp_tbl.frontChoice)) = -1;
ndr = nnz(hddm_tbl.response==-1)/height(hddm_tbl);
lapse_rate = 0.01;



%%
model_name = "nic10"; % or "nic3", use short & memorable model name relevant to data
chains = dir(fullfile("TraceData", model_name));
chains = chains(3:end);
traces = readtable(fullfile(chains(1).folder, chains(1).name));
for i=2:length(chains)
    traces = [traces; readtable(fullfile(chains(i).folder, chains(i).name))];
end
names = traces.Properties.VariableNames;
a_trace = table2array(traces(:, startsWith(names, 'a_subj')));
v_trace = table2array(traces(:, startsWith(names, 'v_subj')));
n_trace = table2array(traces(:, startsWith(names, 't_subj')));
l_trace = table2array(traces(:, startsWith(names, 'alpha_subj')));
b_trace = table2array(traces(:, startsWith(names, 'zt_subj')));
f_trace = table2array(traces(:, startsWith(names, 'forg_subj')));
s_trace = table2array(traces(:, startsWith(names, 'surp_subj')));
st_trace = table2array(traces(:, startsWith(names, 'st')));
a_sex = table2array(traces(:, startsWith(names, 'a_stim_subj')));
v_sex = table2array(traces(:, startsWith(names, 'v_stim_subj')));
n_sex = table2array(traces(:, startsWith(names, 't_stim_subj')));
l_sex = table2array(traces(:, startsWith(names, 'alpha_stim_subj')));
b_sex = table2array(traces(:, startsWith(names, 'zt_stim_subj')));
f_sex = table2array(traces(:, startsWith(names, 'forg_stim_subj')));
s_sex = table2array(traces(:, startsWith(names, 'surp_stim_subj')));
z_trace = table2array(traces(:, strcmp(names, 'z')));
qS_trace = table2array(traces(:, startsWith(names, 'qS')));
qL_trace = table2array(traces(:, startsWith(names, 'qL')));

%% Gelman-Rubin
L = height(traces) / length(chains);
for i=2:length(names)
    chain_means = zeros(size(chains));
    within_var = zeros(size(chains));
    for j=1:length(chain_means)
        chain_means(j) = mean(traces{1+L*(j-1):L*j,i});
        within_var(j) = std(traces{1+L*(j-1):L*j,i});
    end
    grand_mean = mean(chain_means);
    B = L/(length(chains)-1)*sum((chain_means-grand_mean).^2);
    W = mean(within_var);
    disp(strcat(names(i), ": ", num2str(((L-1)/L*W+B/L)/W)))
end

%%
sex_params=[traces.a_1_-traces.a_0_,traces.v_1_-traces.v_0_,traces.zt_1_-traces.zt_0_, traces.t_1_-traces.t_0_,traces.alpha_1_-traces.alpha_0_,traces.forg_1_-traces.forg_0_,traces.surp_1_-traces.surp_0_];
sex_stds=[traces.a_std,traces.v_std,traces.zt_std,traces.t_std,traces.alpha_std,traces.forg_std,traces.surp_std];
labels = ["Boundary Separation","Drift Rate","Bias","Non-Decision Time","Learning Rate", "Forgetfulness", "Surprise"];
effect_size=0.1;
figure
hold on
for i=1:size(sex_params,2)
    normed=sex_params(:,i)./sex_stds(:,i);
    sorted=sort(normed);
    l=sorted(0.025*length(sorted));
    h=sorted(0.975*length(sorted));
    plot_95kd(normed,'c',[247,119,110]/255,@(x, M)0.9*x/M-(i-1))
    disp(labels(i))
    disp(nnz(abs(normed)<effect_size)/length(normed)*100)
    disp(nnz(normed>0)/length(normed)*100)
end
patch('XData',effect_size*[-1,1,1,-1],'YData',[1,1,1-size(sex_params,2),1-size(sex_params,2)],'EdgeColor','none','FaceColor',[0.9,0.9,0.9],'FaceAlpha',0.3)
yticks(1-size(sex_params,2):0)
yticklabels(flip(labels))
ylim([0.5-size(sex_params,2), 1])
xlim([-2,2])
xticks(-2:2)
xlabel('Effect Size')
set(gca,'fontsize',18)
% view([-90,90])

%%
hddm_tbl = readtable('hddm_nic_3sec.csv');
T = height(traces);
N = height(hddm_tbl);
S = length(unique(hddm_tbl.split_by));
subjects = hddm_tbl.subj_idx([1;find(diff(hddm_tbl.split_by)~=0)+1]);

slices = 32;
slice = reshape(1:T, [], slices)';
RT = cell(slices,1);
v_full = cell(slices,1);
a_full = cell(slices,1);
n_tf = zeros(N,1);
conf_tf = zeros(N,1);
choices = cell(slices,1);
conf_full = cell(slices,1);
conf_hat_full = cell(slices,1);
correct_full = cell(slices,1);
rule = cell(slices,1);
rule_num = cell(slices,1);
Vs = cell(slices,1);
Vs_hat = cell(slices,1);

%%
ndr = nnz(hddm_tbl.response==-1)/height(hddm_tbl);
lapse_rate = 0.01;
maxRT = 10;
T = height(traces);
N = height(hddm_tbl);
S = length(unique(hddm_tbl.split_by));
subjects = hddm_tbl.subj_idx([1;find(diff(hddm_tbl.split_by)~=0)+1]);

slices = 32;
slice = reshape(1:T, [], slices)';
RT = cell(slices,1);
v_full = cell(slices,1);
a_full = cell(slices,1);
choices = cell(slices,1);
correct_full = cell(slices,1);
rule = cell(slices,1);
rule_num = cell(slices,1);
Vs = cell(slices,1);
Vs_hat = cell(slices,1);
b_full = cell(slices,1);

rng(626)
%parpool(12) % can get rid of parpool but it helps with efficiency 
for s=1:slices %parfor instead of for if using parpool
    inds = slice(s,:);
    a_temp = a_trace(inds,:);
    l_temp = l_trace(inds,:);
    n_temp = n_trace(inds,:);
    v_temp = v_trace(inds,:);
    b_temp = b_trace(inds,:);
    f_temp = f_trace(inds,:);
    s_temp = s_trace(inds,:);
    st_temp = st_trace(inds);
    z_temp = z_trace(inds);
    qs_temp = qS_trace(inds);
    ql_temp = qL_trace(inds);
    a_stim_temp = a_sex(inds,:);
    l_stim_temp = l_sex(inds,:);
    n_stim_temp = n_sex(inds,:);
    v_stim_temp = v_sex(inds,:);
    b_stim_temp = b_sex(inds,:);
    f_stim_temp = f_sex(inds,:);
    s_stim_temp = s_sex(inds,:);
    Vs_temp = zeros(N,T/slices,2);
    Vs_hat_temp = zeros(N,T/slices,2);
    rn_temp = ones(T/slices,S);
    V = [0,0];
    V_hat = [0,0];
    v_t = zeros(N,T/slices);
    a_t = zeros(N,T/slices);
    b_t = zeros(N,T/slices);
    v_that = zeros(N,T/slices);
    b_that = zeros(N,T/slices);
    n_t = zeros(N,T/slices);
    choices_temp = zeros(T/slices,N);
    rule_temp = false(T/slices,N);
    rt_temp = zeros(T/slices, N);
    correct_temp = false(T/slices, N);
    Q = [qs_temp(1), qs_temp(1), ql_temp(1)];
    Q_hat = [qs_temp(1), qs_temp(1), ql_temp(1)];
    ses = 1;
    blockNum = 1;
    for j=1:T/slices
        for i=1:N
            if i == 1 || hddm_tbl.split_by(i) ~= hddm_tbl.split_by(i-1)
                if i~= 1
                    rn_temp(j,ses) = rn_temp(j,ses) + (blockNum - 1) / 5;
                end
                Q = [qs_temp(j), qs_temp(j), ql_temp(j)];
                Q_hat = [qs_temp(j), qs_temp(j), ql_temp(j)];
                blockNum = 1;
                ses = hddm_tbl.split_by(i);
                rules = 1:3;
                curRule = randsample(length(rules), 1);
                rules(curRule) = [];
            end
            alpha_t = a_temp(j, hddm_tbl.subj_idx(i));
            lr_t = logit(l_temp(j, hddm_tbl.subj_idx(i)));
            delta_t = v_temp(j, hddm_tbl.subj_idx(i));
            ndt_t = n_temp(j, hddm_tbl.subj_idx(i));
            bias_t = b_temp(j, hddm_tbl.subj_idx(i));
            forg_t = logit(f_temp(j, hddm_tbl.subj_idx(i)));
            surp_t = exp(s_temp(j, hddm_tbl.subj_idx(i)));
            n_t(i,j) = ndt_t;
            a_t(i,j) = alpha_t;
            if hddm_tbl.rt(i) ~= -1
                V(1) = Q(1) + (hddm_tbl.light(i) == 1) * (Q(3));
                V(2) = Q(2) + (hddm_tbl.light(i) == 2) * (Q(3));
                Vs_temp(i,j,:) = V;
                v_t(i,j) = (V(1) - V(2)) * delta_t;
                b_t(i,j) = logit((Q(1) - Q(2)) * bias_t + z_temp(j));
                pe = hddm_tbl.feedback(i) - Q(hddm_tbl.response(i));
                adj_lr = lr_t * abs(pe) ^ surp_t;
                Q(hddm_tbl.response(i)) = Q(hddm_tbl.response(i)) + adj_lr * pe;
                Q(abs(hddm_tbl.response(i)-3)) = Q(abs(hddm_tbl.response(i)-3)) * (1-forg_t);
                if hddm_tbl.response(i) == hddm_tbl.light(i)
                    pe = hddm_tbl.feedback(i) - Q(3);
                    adj_lr = lr_t * abs(pe) ^ surp_t;
                    Q(3) = Q(3) + adj_lr * pe;
                else
                    Q(3) = Q(3) * (1-forg_t);
                end
            end
            if rand < ndr
                choices_temp(j,i) = -1;
                rt_temp(j,i) = -1;
                correct_temp(j,i) = 0;
            else
                V_hat(1) = Q_hat(1) + (hddm_tbl.light(i) == 1) * (Q_hat(3));
                V_hat(2) = Q_hat(2) + (hddm_tbl.light(i) == 2) * (Q_hat(3));
                Vs_hat_temp(i,j,:) = V_hat;
                v_that(i,j) = (V_hat(1) - V_hat(2)) * delta_t;
                b_that(i,j) = logit((Q_hat(1) - Q_hat(2)) * bias_t + z_temp(j));
                if rand < lapse_rate
                    choices_temp(j,i) = -(round(rand)-2);
                    rt_temp(j,i) = unifrnd(0,maxRT);
                else
                    rt_temp(j,i) = wienerrng(a_t(i,j), n_t(i) + st_temp(j) * (rand - 0.5),b_that(i,j)*a_t(i,j),v_that(i,j));
                    choices_temp(j, i) = -((rt_temp(j,i)>0)-2);
                    rt_temp(j,i) = abs(rt_temp(j,i));
                end
                if choices_temp(j,i) == curRule || (choices_temp(j,i) == hddm_tbl.light(i) && 3 == curRule)
                    correct_temp(j,i) = 1;
                else
                    correct_temp(j,i) = 0;
                end
                pe = correct_temp(j,i) - Q_hat(choices_temp(j,i));
                adj_lr = lr_t * abs(pe) ^ surp_t;
                Q_hat(choices_temp(j,i)) = Q_hat(choices_temp(j,i)) + adj_lr * pe;
                Q_hat(abs(choices_temp(j,i)-3)) = Q_hat(abs(choices_temp(j,i)-3)) * (1-forg_t);
                if choices_temp(j,i) == hddm_tbl.light(i)
                    pe = correct_temp(j,i) - Q_hat(3);
                    adj_lr = lr_t * abs(pe) ^ surp_t;
                    Q_hat(3) = Q_hat(3) + adj_lr * pe;
                else
                    Q_hat(3) = Q_hat(3) * (1-forg_t);
                end
            end
            rule_temp(j,i) = curRule > 2;
            if correct_temp(j,i)
                blockNum = blockNum + 1;
                if blockNum == 6
                    if isempty(rules)
                        rules = 1:3;
                    end
                    curRule = rules(randsample(length(rules), 1));
                    rules(rules==curRule) = [];
                    rn_temp(j,ses) = rn_temp(j,ses) + 1;
                end
            else
                blockNum = 1;
            end
            blockNum = mod(blockNum - 1, 5) + 1;
        end
        rn_temp(j,end) = rn_temp(j,end) + (blockNum - 1) / 5;
        disp(j)
    end
    RT{s} = rt_temp;
    rule_num{s} = rn_temp;
    choices{s} = choices_temp;
    correct_full{s} = correct_temp;
    v_full{s} = v_t;
    a_full{s} = a_t;
    b_full{s} = b_t;
    Vs{s} = Vs_temp;
    Vs_hat{s} = Vs_hat_temp;
    rule{s} = rule_temp;
end

delete(gcp('nocreate'))
b_full_mat = cell2mat(b_full');
RT_mat=cell2mat(cellfun(@(x) transpose(x),RT,'UniformOutput',false)');
dV_full_mat = cell2mat(cellfun(@(x) x(:,:,1)',Vs,'UniformOutput',false)) - cell2mat(cellfun(@(x) x(:,:,2)',Vs,'UniformOutput',false));
dV_hat_full_mat = cell2mat(cellfun(@(x) x(:,:,1)',Vs_hat,'UniformOutput',false)) - cell2mat(cellfun(@(x) x(:,:,2)',Vs_hat,'UniformOutput',false));
choices_mat = cell2mat(cellfun(@(x) transpose(x),choices,'UniformOutput',false)');
rule_num_mat = cell2mat(rule_num);
correct_full_mat = cell2mat(cellfun(@(x) transpose(x),correct_full,'UniformOutput',false)');
rule_mat = cell2mat(cellfun(@(x) transpose(x),rule,'UniformOutput',false)');

%%
figure
subplot(2,2,1)
points = -0.75:0.001:0.75;
choice_curves = zeros(length(points),T);
choice_curves_hat = zeros(length(points),T);
rng(626)
for i=1:T
    data_choice=hddm_tbl.response(hddm_tbl.rt>0);
    [temp,I]=sort(dV_full_mat(i,hddm_tbl.rt>0)+rand(1,nnz(hddm_tbl.rt>0))*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(data_choice(I),0.1,'samplepoints',temp),points);
    hat_rt=choices_mat(RT_mat(:,i)>0,i);
    [temp,I]=sort(dV_hat_full_mat(i,RT_mat(:,i)>0)+rand(1,nnz(RT_mat(:,i)>0))*10e-6);
    choice_curves_hat(:,i)=interp1(temp,movmean(hat_rt(I),0.1,'samplepoints',temp),points);
    if mod(i,100)==0
        disp(i)
    end
end
mcurve = 100*nanmedian(-choice_curves+2,2);
sor_curve = 100*(sort(-choice_curves+2,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
mcurve2=100*nanmedian(-choice_curves_hat+2,2);
sor_curve2 = 100*(sort(-choice_curves_hat+2,2));
lcurve2 = sor_curve2(:,0.025*4000);
hcurve2 = sor_curve2(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
patch('XData',[points,flip(points)], 'YData', [lcurve2;flip(hcurve2)],'FaceColor',[255,154,0]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'-','Color',[180,0,255]/255,'LineWidth',2)
plot(points,mcurve2,'-','Color',[255,154,0]/255,'LineWidth',2)
xlabel('\DeltaV')
ylabel('% Left Choice')
xlim([points(1),points(end)])
xticks(-1:0.5:1)
yticks(0:50:100)
ylim([0,100])
set(gca,'fontsize',18)
%%
subplot(2,2,2)
points = -0.75:0.001:0.75;
choice_curves = zeros(length(points),T);
choice_curves_hat = zeros(length(points),T);
rng(626)
for i=1:T
    data_rt=hddm_tbl.rt(hddm_tbl.rt>0);
    [temp,I]=sort(dV_full_mat(i,hddm_tbl.rt>0)+rand(1,nnz(hddm_tbl.rt>0))*10e-6);
    choice_curves(:,i)=interp1(temp,movmean(data_rt(I),0.05,'samplepoints',temp),points);
    hat_rt=RT_mat(RT_mat(:,i)>0,i);
    [temp,I]=sort(dV_hat_full_mat(i,RT_mat(:,i)>0)+rand(1,nnz(RT_mat(:,i)>0))*10e-6);
    choice_curves_hat(:,i)=interp1(temp,movmean(hat_rt(I),0.05,'samplepoints',temp),points);
    if mod(i,100)==0
        disp(i)
    end
end
mcurve = nanmedian(choice_curves,2);
sor_curve = (sort(choice_curves,2));
lcurve = sor_curve(:,0.025*4000);
hcurve = sor_curve(:, 0.975*4000);
mcurve2=nanmedian(choice_curves_hat,2);
sor_curve2 = (sort(choice_curves_hat,2));
lcurve2 = sor_curve2(:,0.025*4000);
hcurve2 = sor_curve2(:, 0.975*4000);
hold on
patch('XData',[points,flip(points)], 'YData', [lcurve;flip(hcurve)],'FaceColor',[180,0,255]/255,'EdgeColor','none','FaceAlpha',0.3)
patch('XData',[points,flip(points)], 'YData', [lcurve2;flip(hcurve2)],'FaceColor',[255,154,0]/255,'EdgeColor','none','FaceAlpha',0.3)
plot(points,mcurve,'-','Color',[180,0,255]/255,'LineWidth',2)
plot(points,mcurve2,'-','Color',[255,154,0]/255,'LineWidth',2)
xlabel('\DeltaV')
ylabel('RT (s)')
xlim([points(1),points(end)])
xticks(-1:0.5:1)
%yticks(0.55:0.1:0.75)
%ylim([0.55,0.75])
set(gca,'fontsize',18)
%%
subplot(2,2,3)
hold on
pts = 0:0.001:10;
all_vals = zeros(size(pts));
bw = 0.05;
for s = 1:slices
    for i=1:50:T/slices
        temp = RT{s}(i,:);
        temp = squeeze(temp);
        temp=rmmissing(temp(temp>0));
        if ~isempty(temp)
            vals = ksdensity(temp,pts,'bandwidth',bw);
            plot(pts, vals, "Color",[0.949,0.631,0.008,0.05])
            all_vals = all_vals + vals;
        end
    end
end
hold on
vals=ksdensity(squeeze(hddm_tbl.rt(hddm_tbl.rt>-1)),pts,'bandwidth',bw);
plot(pts, vals, "Color",[180,0,255]/255,'LineWidth',1.5)
xticks(0:10)
xlim([0,10])
xlabel('RT (s)')
yticks([])
set(gca,'fontsize',18)
set(gca,'linewidth',2)
set(gca,'fontname','Helvetica')

%%
temp_tbl = readtable('temp_tbl.csv');
subplot(2,2,4)
sessions = unique(hddm_tbl.split_by);
acc_side = zeros(length(sessions), 4000);
acc_light = zeros(length(sessions), 4000);
rat_side = zeros(size(sessions));
rat_light = zeros(size(sessions));
for i=1:length(sessions)
    rule_seg = rule_mat(hddm_tbl.split_by==i,:);
    correct_seg = correct_full_mat(hddm_tbl.split_by==i,:);
    acc_side(i,:) = sum(correct_seg.*~rule_seg)./sum(~rule_seg);
    acc_light(i,:) = sum(correct_seg.*rule_seg)./sum(rule_seg);
    tbl_seg = temp_tbl(hddm_tbl.split_by==i,:);
    rat_side(i) = mean(tbl_seg.performance(~tbl_seg.LightRule));
    rat_light(i) = mean(tbl_seg.performance(tbl_seg.LightRule==1));
end
hold on
pts = 0:0.005:1;
for i=1:50:4000
    plot(pts,ksdensity(acc_side(:,i),pts), 'Color',[0.31,0.73,0.90,0.05])
end
for i=1:50:4000
    plot(pts,ksdensity(acc_light(:,i),pts), 'Color', [0.83,0.07,0.82,0.05])
end
plot(pts,ksdensity(rat_side, pts),'Color',[0.31,0.73,0.90],'LineWidth',2)
plot(pts,ksdensity(rat_light, pts),'Color',[0.83,0.07,0.82],'LineWidth',2)
set(gca,'fontsize',18)
xlabel('Accuracy')
yticks([])
ylim([0,8])
set(gca,'linewidth',2)
set(gca,'fontname','Helvetica')

%% this part still has errors but I have better ways to visualize rules anyway
figure
sessions = unique(temp_tbl.Session);
rules_completed = zeros(size(sessions));
for i=1:length(sessions)
    subtbl = temp_tbl(strcmp(temp_tbl.Session, sessions(i)), :);
    rules_completed(i) = length(find(~strcmp(subtbl.rule(2:end), subtbl.rule(1:end-1))));
    rules_completed(i) = rules_completed(i) + (subtbl.block(end) - 1) / 5;
end
plot_95kd(mean(rule_num_mat,2)-1, 'c',[247,119,110]/255)
yticks([])
xline(mean(rules_completed),'k--','LineWidth',2)
xlabel('Rules Completed')
set(gca,'fontsize',18)
set(gca,'linewidth',2)

%%
function plot_95kd(data,c1,c2,warp)
    [ad, ax] = ksdensity(data,'NumPoints',1000);
    sorted=sort(data);
    l=sorted(0.025*length(sorted));
    h=sorted(0.975*length(sorted));
    if nargin > 3
        patch('XData',[ax(1),ax,ax(end)],'YData',warp([0,ad,0], max(ad)),'EdgeColor','none','FaceColor',c1)
        patch('XData',[l,ax(ax>l&ax<h),h],'YData',warp([0,ad(ax>l&ax<h),0], max(ad)),'EdgeColor','none','FaceColor',c2)
        [~,I] = min(abs(ax-median(data)));
        plot(median(data)*[1,1],warp(ad(I)*[0,1],max(ad)),'Color',0.8*c2,'linewidth',2)
    else
        patch('XData',[ax(1),ax,ax(end)],'YData',[0,ad,0],'EdgeColor','none','FaceColor',c1)
        patch('XData',[l,ax(ax>l&ax<h),h],'YData',[0,ad(ax>l&ax<h),0],'EdgeColor','none','FaceColor',c2)
    end
end