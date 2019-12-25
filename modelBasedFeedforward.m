function [feedforwardController,forwardOrder] = modelBasedFeedforward(G,ffMethod)
[~,~,Ts] = tfdata(G,'v');

ffMethod = lower(ffMethod);

switch ffMethod
    case 'zpetc'
        feedforwardController = ZPETC(G);
    case 'zmetc'
        feedforwardController = ZMETC(G);
    case 'ignore'
        feedforwardController = nonminimumIgnore(G);
    case 'seriestruncation'
    otherwise
        
        
        
end

forwardOrder = numel(zero(feedforwardController)) - numel(pole(feedforwardController));


    function F = ZPETC(G)
        %%
        [z,p,k,Ts] = zpkdata(G,'v');
        % relativeDegree = numel(p) - numel(z);
        bound = 1;
        index = abs(z) >= bound;
        nonminimumZero = z( index );
        minimumZero = z( ~index );
        F = zpk(p,minimumZero,1/k,Ts);
        %%
        %         for i = 1:numel(nonminimumZero)
        %             temp = tf( [-1 * nonminimumZero(i),1] ,1 ,Ts,'variable','z^-1');
        %             F = F * temp * (1 / (abs(1 - nonminimumZero(i))).^2 );
        %         end
        i = 1;
        while i <= numel(nonminimumZero)
            if(isreal(nonminimumZero(i)))
                temp = tf( [-1 * nonminimumZero(i),1] ,1 ,Ts,'variable','z^-1');
                F = F * temp * (1 / (abs(1 - nonminimumZero(i))).^2 );
                i = i + 1;
            else
                temp = tf( [ ( abs(nonminimumZero(i)) ).^2 ,-2 * real(nonminimumZero(i)) ,1] ,1 ,Ts,'variable','z^-1');
                F = F * temp * (1 / ( 1 + ( abs(nonminimumZero(i)) ).^2 -2 * real(nonminimumZero(i))).^2);
                i = i + 2;
            end
            
        end
        [b,a,T] = tfdata(F,'v');
        F = tf(b,a,T,'variable',F.Variable);
        
    end

    function F = ZMETC(G)
        %%
        [z,p,k,Ts] = zpkdata(G,'v');
        % relativeDegree = numel(p) - numel(z);
        bound = 1;
        index = abs(z) >= bound;
        nonminimumZero = z( index );
        minimumZero = z( ~index );
        F = zpk(p,minimumZero,1/k,Ts);
        %%
        %         for i = 1:numel(nonminimumZero)
        %             temp = tf(1 , [-1 * nonminimumZero(i),1] ,Ts,'variable','z^-1');
        %             F = F * temp;
        %         end
        
        i = 1;
        while i <= numel(nonminimumZero)
            if(isreal(nonminimumZero(i)))
                temp = tf( 1 ,[-1 * nonminimumZero(i),1] ,Ts,'variable','z^-1');
                F = F * temp;
                i = i + 1;
            else
                temp = tf( 1 ,[ ( abs(nonminimumZero(i)) ).^2 ,-2 * real(nonminimumZero(i)) ,1] ,Ts,'variable','z^-1');
                F = F * temp ;
                i = i + 2;
            end
            
        end
        z = tf('z',Ts);
        F = minreal(F * z^( -1 * 2 * numel(nonminimumZero)));
        [b,a,T] = tfdata(F,'v');
        F = tf(b,a,T,'variable',F.Variable);
    end

    function F = nonminimumIgnore(G)
        %%
        [z,p,k,Ts] = zpkdata(G,'v');
        % relativeDegree = numel(p) - numel(z);
        bound = 0.99999999999;
        index = abs(z) > bound;
        nonminimumZero = z( index );
        minimumZero = z( ~index );
        F = zpk(p,minimumZero,1/k,Ts);
        for i = 1:numel(nonminimumZero)
            %             temp = tf( [-1 * nonminimumZero(i),1] ,1 ,Ts,'variable','z^-1');
            F = F * (1 / (abs(1 - nonminimumZero(i))) );
        end
        z = tf('z',Ts);
        F = minreal(F * z^( -1 *  numel(nonminimumZero)));
    end

    function F = seriesTruncation(G)
        
    end


end