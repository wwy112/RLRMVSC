function [result] = CluEst(Y,C,truth)
% Y -- low-dimensional representation truth-- true label
Maxm = 10;
ACCi = zeros(Maxm,1);
NMIi = zeros(Maxm,1);
fscore = zeros(Maxm,1);
precision = zeros(Maxm,1);
recall = zeros(Maxm,1);
ARI = zeros(Maxm,1);
RI = zeros(Maxm,1);
purity = zeros(Maxm,1);
for m = 1 : Maxm
    grps = kmeans(Y,C,'EmptyAction','drop');
    P_label = bestMap(truth, grps);
    ACCi(m) = length(find(truth == P_label))/length(truth);
    [~, NMIi(m), ~] = compute_nmi(truth, grps);
    [fscore(m), precision(m), recall(m)] = compute_f(truth,grps);
    [ARI(m),RI(m),~,~] = RandIndex(truth,grps);
    purity(m) = pur_fun(truth,P_label);
end
result = [mean(ACCi),std(ACCi);mean(NMIi),std(NMIi);mean(precision),std(precision);mean(recall),std(recall);mean(fscore),std(fscore);mean(purity),std(purity);mean(ARI),std(ARI);mean(RI),std(RI)];
