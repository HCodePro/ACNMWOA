%% 反向学习
function [ NewPop, FEs, NewFit ] = GOL( Pop,Fit,dim,ub,lb,thdim,fname,Pxy,FEs,type)
%   Pop种群，N种群大小，dim问题维度，lb下界，ub上界，fobj适应度函数，type反向学习模式选择参数，0为固定边界，1位动态边界。
%   默认上下边界向量是行方向，如果实际情况是列方向，请将上下边界向量转置后使用
	if type==0 %%一般随机反向学习
		for i = 1:size(Pop,1)
			 for j = 1:size(Pop,2)
				TemPop(i,j) = rand()*(ub(j)+lb(j)) - Pop(i,j);
				if TemPop(i,j) < lb(j) || TemPop(i,j) > ub(j)
					TemPop(i,j) = lb(j) + rand() * (ub(j) - lb(j));  
				end
			 end
			 [TemFit(1,i),FEs,TemPop(i,:)] =sigle_evaluation(TemPop(i,:),dim,thdim,fname,Pxy,FEs);
		end
		ComPop = [Pop; TemPop];
		ComFit = [Fit, TemFit];
		[~,index] = sort(ComFit,'descend');
		for i = 1:size(Pop,1)
			NewPop(i,:) = ComPop(index(i),:);
			NewFit(1,i) = ComFit(1, index(i));
		end
	elseif type==1 %%动态边界随机反向学习
		dub=max(Pop);
		dlb=min(Pop);
		for i = 1:size(Pop,1)
			 for j = 1:size(Pop,2)
				TemPop(i,j) = rand()*(dub(j)+dlb(j)) - Pop(i,j);
				if TemPop(i,j) < dlb(j) || TemPop(i,j) > dub(j)
					TemPop(i,j) = dlb(j) + rand() * (dub(j) - dlb(j));  
				end
			 end
			 [TemFit(1,i),FEs,TemPop(i,:)] =sigle_evaluation(TemPop(i,:),dim,thdim,fname,Pxy,FEs);
		end
		ComPop = [Pop; TemPop];
		ComFit = [Fit, TemFit];
		[~,index] = sort(ComFit,'descend');
		for i = 1:size(Pop,1)
			NewPop(i,:) = ComPop(index(i),:);
			NewFit(1,i) = ComFit(1, index(i));
		end		
	end
end

