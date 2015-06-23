phiRes = 5;
phiDotRes = 5;
% M = zeros(phiRes*phiDotRes,1);
for i = 1:phiRes
	for j = 1:phiDotRes
		load(['ensemble_traj_10_',num2str(i),'_',num2str(j),'.mat']);
		figErrBar = errorbar(HVec,avgChi0Vec,stdChi0Vec,'-ob');
        title(['Initial condition:(',num2str(xInit(1)),...
            ',',num2str(xInit(2)),')'])
		hold on
        axis([1 11 0 2])
        pause(1)
        M((i-1)*phiDotRes + j) = getframe;
	end
end
