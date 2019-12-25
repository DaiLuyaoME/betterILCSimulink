function [compensationPart,forwardOrder,coef] = extraCompensation(G,ffMethod)
[~,~,Ts] = tfdata(G,'v');

ffMethod = lower(ffMethod);

switch ffMethod
    case 'zpetc'
        [compensationPart,coef] = ZPETC(G);
    case 'zmetc'
        [compensationPart,coef] = ZMETC(G);
    case 'ignore'
        [compensationPart,coef] = nonminimumIgnore(G);
    case 'seriestruncation'
    otherwise
        
        
        
end

forwardOrder = numel(zero(compensationPart)) - numel(pole(compensationPart));


    function [F,coef] = ZPETC(G)
        %%
        [z,p,k,Ts] = zpkdata(G,'v');
        % relativeDegree = numel(p) - numel(z);
        bound = 1;
        index = abs(z) >= bound;
        nonminimumZero = z( index );
        alpha1 = sum( nonminimumZero./(1 - nonminimumZero).^2 );
        z = tf('z',Ts);
        F = z * ((1-z^-1)/Ts)^2;
        coef = alpha1 * Ts;
        
        %
        %
        %         [b,a,T] = tfdata(F,'v');
        %         F = tf(b,a,T,'variable',F.Variable);
        
    end

    function [F,coef] = ZMETC(G)
        %%
        [z,p,k,Ts] = zpkdata(G,'v');
        % relativeDegree = numel(p) - numel(z);
        bound = 1;
        index = abs(z) >= bound;
        nonminimumZero = z( index );
        alpha1 = -1 * sum( 2 * nonminimumZero ./ (1 - nonminimumZero) );
        z = tf('z',Ts);
        F =  z * (1-z^-1)/Ts;
        coef = alpha1 * Ts;
        
        
        
        %         [b,a,T] = tfdata(F,'v');
        %         F = tf(b,a,T,'variable',F.Variable);
    end

    function [F,coef] = nonminimumIgnore(G)
        %%
        [z,p,k,Ts] = zpkdata(G,'v');
        % relativeDegree = numel(p) - numel(z);
        bound = 1;
        index = abs(z) >= bound;
        nonminimumZero = z( index );
        alpha1 = -1 * sum( nonminimumZero ./ (1 - nonminimumZero) );
        z = tf('z',Ts);
        F = (1-z^-1)/Ts;
        coef = alpha1 * Ts;
    end

    function F = seriesTruncation(G)
        
    end


end