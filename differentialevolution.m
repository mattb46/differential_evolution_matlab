function [xOut,fOut,x,fVal]  = differentialevolution(funEval, xL, xU, Np, gMax, Cr, F)
  % description
  % 1) generate random test vectors within bounds
  % 2) generate trial vectors
  %     a) add a weighted difference of two random test vectors with another
  %        random test vector
  %     b) crossover to enhance diversity, replace randomly selected indices
  %        of the trial vectors with the test vectors
  % 3) compare test and trial vectors 1:1 and replace inferior test vectors
  %    for the next generation
  % 4) repeat 2-3
  %
  % variables
  % x     ->  target vector
  % u     ->  trial vector
  % xL    ->  lower bound
  % xU    ->  upper bounds
  % Np    ->  population vector number
  % D     ->  number of dimensions
  % gMax  ->  max # of generations
  % fVal  ->  function evaluation of test vector
  % F     ->  difference vector scale/weighted
  % Cr    ->  between 0-1, 0 replacing 1 vector, 1 replacing
  %           all vectors
  %
  D = length(xL);
  x = ones(D,Np,gMax);
  u = ones(D,Np,gMax);
  fVal = ones(Np,gMax);
  fValU = ones(Np,gMax);
  g=1;
  for i=1:Np % initialize population
    for j=1:D
      x(j,i,g) = xL(j)+unifrnd(0,1)*(xU(j)-xL(j));
    end
    fVal(i,1) = feval(funEval,x(:,i,g)); % Evaluate and store f(x(i),g)
  end
  for g=1:gMax % generation loop
    for i=1:Np % generate a trial population
      jrandi = floor(unifrnd(0,1)*D)+1; % randomly selet a parameter
      r1=floor(unifrnd(0,1)*Np)+1;
      while (r1==i)
        r1=floor(unifrnd(0,1)*Np)+1;
      end
      r2=floor(unifrnd(0,1)*Np)+1;
      while(r2==i || r2==r1)
        r2=floor(unifrnd(0,1)*Np)+1;
      end
      r3=floor(unifrnd(0,1)*Np)+1;
      while (r3==i || r3==r2 || r3==r1)
        r3=floor(unifrnd(0,1)*Np)+1;
      end
      for j=1:D % generate a trial vector
        if(unifrnd(0,1)<=Cr || j==jrandi)
          u(j,i,g) = x(j,r1,g)+F*(x(j,r2,g)-x(j,r3,g));
        else
          u(j,i,g) = x(j,i,g);
        end
      end
    end
    for i=1:Np
      fValU(i,g) = feval(funEval,u(:,i,g));
    end
    for i=1:Np % select new population
      if (fValU(i,g) <= fVal(i,g)) % evaluate trial vector and compare w/target vector % f(u(i,g) <= f(x(i,g)
        for j=1:D
          x(j,i,g+1) = u(j,i,g); % replace inferior target
        end
        fVal(i,g+1) = fValU(i,g);
      else
          for j=1:D
            x(j,i,g+1) = x(j,i,g);
          end
          fVal(i,g+1) = fVal(i,g);
      end
    end
  end
  fOut = fVal(i,gMax);
  xOut = x(:,i,gMax);
  for i=2:Np
    if(fVal(i,gMax)<fOut)
      fOut = fVal(i,gMax);
      xOut = x(:,i,gMax);
    end
  end
end
