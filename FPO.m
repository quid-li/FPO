function [bestF, bestX, Curve] = FPO(N, T, lb, ub, dim, fobj) 
	X = lb + (ub - lb) .* rand(N, dim); 
	newX = zeros(N, dim); 
    Fit = zeros(N, 1); 
	bestX = X(end, :);
    bestF = inf;
    for i = 1:N
        Fit(i) = fobj(X(i, :)); 
		if Fit(i) < bestF
			bestX = X(i, :);
			bestF = Fit(i);
		end
    end
	ratios = [3, 2, 1, 3, 3];
	seasons = [0, cumsum(floor(T * ratios(1:4) / 12)), T];
	season_index = 1; 
	RM = [8, 2, 1, 3, 5; 5, 8, 2, 1, 3; 3, 5, 8, 2, 1; 1, 3, 5, 8, 2; 2, 1, 3, 5, 8];
	subpop_sizes = ceil(N * RM(5, :) / sum(RM(5, :)));
	subpop_sizes(5) = subpop_sizes(5) - (sum(subpop_sizes) - N); 
	GR = 2.0/(1 + sqrt(5));
    CRmax = 0.99;
    CRmin = 0.4;
	Curve = zeros(1, T); 
    for t = 0:T-1
		if t == seasons(season_index)
			subpop_sizes = circshift(subpop_sizes, 1);  
			cs = cumsum(subpop_sizes);
			season_index = season_index + 1;
        end
        P = t / T;
		CR = CRmin + (CRmax - CRmin) * P;
		DP = (1 - P)^2; 
		F = 0.1 + 0.5 * DP;
        RT = F * trnd(1, N, dim); 
		RU = rand(N, dim); 
		[~, sorted_idx] = sort(Fit, 'descend');
		X = X(sorted_idx, :);
		Fit = Fit(sorted_idx);
		worstX = X(1, :); 
		centerX = mean(X, 1);
		for i = 1 : N
			R = randi([1, N], 1, 2);
			if i <= cs(1)
				newX(i, :) = X(i, :) + RU(i, :) .* (X(R(1), :) - randi(2) * X(i, :));
			elseif i <= cs(2)
				newX(i, :) = bestX + GR * (bestX - X(i, :)) + DP * RU(i, :) .* (X(R(2), :) - worstX);
			elseif i <= cs(3)
				mask = RU(i, :) < P;
				newX(i, :) = mask .* bestX + (~mask) .* X(i, :) + RT(i, :) .* (centerX - X(i, :));
			elseif i <= cs(4)
				mask = RU(i, :) < CR;
				newX(i, :) = mask .* bestX + (~mask) .* X(i, :) + RT(i, :) .* (X(R(1), :) - X(R(2), :));
			else
				mask = RU(i, :) < CR;   
				newX(i, :) = X(i, :) + mask .* (F * (bestX - X(i, :)) + RT(i, :) .* (X(R(1), :) - X(R(2), :)));
			end
			violate = (newX(i, :) - lb < 0) | (ub - newX(i, :) < 0);
			if any(violate)
				mutant = X(R(2), :) + RU(i, :) .* (X(R(1), :) - X(R(2), :));
				newX(i, violate) = mutant(violate);
			end
			newF = fobj(newX(i, :));
			if newF < Fit(i)
				X(i, :) = newX(i, :);
				Fit(i) = newF;
				if newF < bestF
					bestX = newX(i, :);
					bestF = newF;
				end
			end
		end
		Curve(t+1) = bestF;                      
    end
end
