
function gradient_st(tree::Matrix{Float64},x::Matrix{Float64},var::Int64)::Vector{Float64}
    N::Int64 = size(x)[1]
    terminalnodes::Vector{Int64} = find(x->(x==1),tree[:,6])
    terminal::Matrix{Float64} = tree[terminalnodes,:]
    Nt::Int64 = length(terminalnodes)
    logimat::Matrix{Float64} = Matrix{Float64}(N,Nt)
    node::Vector{Float64} = Vector{Float64}(9)
    logit::Vector{Float64} = Vector{Float64}(N)
    dlogit::Vector{Float64} = Vector{Float64}(N)
    x_split::Vector{Float64} = Vector{Float64}(N)
    variable::Int64 = 1::Int64
    parent::Int64 = 1::Int64
    logitaux::Vector{Float64}=Vector{Float64}(N)
    dlogitaux::Vector{Float64}=Vector{Float64}(N)
    for i in 1:Nt
        node=terminal[i,:] 
        variable = floor(Int,node[7])
        x_split = x[:,variable]
        [logit[j] = 1/(1+exp(-node[4]*(x_split[j]-node[3]))) for j in 1:N]
        if variable == var
            dlogit = node[4]*logit.*(1-logit)
        else
            dlogit = zeros(N)
        end
        if node[1] == 2.0
            logit = 1.0-logit
            dlogit = -dlogit
        end

        parent = floor(Int,node[5])
        while parent!=0
            node = tree[parent,:]
            variable = floor(Int,node[7])
            x_split = x[:,variable]
            [logitaux[j] = 1/(1+exp(-node[4]*(x_split[j]-node[3]))) for j in 1:N]
            if variable == var
                dlogitaux = node[4]*logitaux.*(1-logitaux)
            else
                dlogitaux = zeros(N)
            end
            if node[1] == 2.0
                logitaux = 1.0-logitaux
                dlogitaux = -dlogitaux
            end
            [dlogit[j] = logit[j]*dlogitaux[j] + dlogit[j]*logitaux[j] for j in 1:N]
            [logit[j] = logit[j]*logitaux[j] for j in 1:N]
            parent = floor(Int,node[5])
        end
        logimat[:,i] = dlogit
    end
    fitted::Vector{Float64} = logimat*terminal[:,2]
    return fitted
end

function fastlm(x::Matrix{Float64},y::Vector{Float64})::Vector{Float64}
    K::Int64 = size(x,2)
    N::Int64 = size(x,1)
    C::Matrix{Float64} = Matrix{Float64}(K,K)
    D::Matrix{Float64} = Matrix{Float64}(K,N)
    b::Vector{Float64} = Vector{Float64}(K)
    ssr::Float64 = 0
    res::Vector{Float64} = Vector{Float64}(K+1)
    fit::Vector{Float64} = Vector{Float64}(N)
    gemm!('T','N',1.0,x,x, 0.0, C)
    gemm!('N','T',1.0,inv(C),x, 0.0, D)
    gemv!('N',1.0,D,y, 0.0, b)
    gemv!('N',1.0,x,b,0.0,fit)
    ssr = sum( [(y[i]-fit[i])^2 for i in 1:N])
    res[1] = ssr
    res[2:end] = b
    return res
end


function eval_tree(x::Matrix{Float64},tree::Matrix{Float64})::Vector{Float64}
    N::Int64 = size(x)[1]
    terminalnodes::Vector{Int64} = find(x->(x==1),tree[:,6])
    terminal::Matrix{Float64} = tree[terminalnodes,:]
    Nt::Int64 = length(terminalnodes)
    logimat::Matrix{Float64} = Matrix{Float64}(N,Nt)
    node::Vector{Float64} = Vector{Float64}(9)
    logit::Vector{Float64} = Vector{Float64}(N)
    x_split::Vector{Float64} = Vector{Float64}(N)
    variable::Int64 = 1::Int64
    parent::Int64 = 1::Int64
    logitaux::Vector{Float64}=Vector{Float64}(N)
    for i in 1:Nt
        node=terminal[i,:] 
        variable = floor(Int,node[7])
        x_split = x[:,variable]
        [logit[j] = 1/(1+exp(-node[4]*(x_split[j]-node[3]))) for j in 1:N]
        if node[1] == 2.0
            logit = 1.0-logit
        end
        parent = floor(Int,node[5])
        
        while parent!=0
            node = tree[parent,:]
            variable = floor(Int,node[7])
            x_split = x[:,variable]
            [logitaux[j] = 1/(1+exp(-node[4]*(x_split[j]-node[3]))) for j in 1:N]
            if node[1] == 2.0
                logitaux = 1.0-logitaux
            end
            [logit[j] = logit[j]*logitaux[j] for j in 1:N]
            parent = floor(Int,node[5])
        end
        logimat[:,i] = logit
    end
    fitted::Vector{Float64} = logimat*terminal[:,2]

    return fitted
end


function grow_tree(x::Matrix{Float64},y::Vector{Float64},p::Float64,d_max::Int64,gamma::Vector{Float64},node_obs::Int64)::Tuple{Array{Float64,2},Array{Float64,1}}
    
    bf::Int64 = 0
    K::Int64 = size(x,2)
    N::Int64 = length(y)
    (tree::Matrix{Float64}, pmat::Matrix{Float64}) = grow0(x,y,p,N,K,gamma,node_obs)::Tuple{Array{Float64,2},Array{Float64,2}}   
    iter::Int64 = 1
    while iter <= d_max
       (tree, pmat, bf, iter) = grow(x,y,p,N,K,gamma,node_obs,tree,pmat,d_max,iter,bf)::Tuple{Array{Float64,2},Array{Float64,2},Int64,Int64} 
    end
   
    middlenodes::Vector{Int64} = find(x->(isnan(x)),pmat[1,:])
    [pmat[i,j] = 0 for i in 1:N, j in middlenodes]
    fitted::Vector{Float64} = pmat*tree[:,2]
    
    return tree, fitted
end



function grow(x::Matrix{Float64},y::Vector{Float64},p::Float64,N::Int64,K::Int64,gamma::Vector{Float64},node_obs::Int64,tree::Matrix{Float64},pmat::Matrix{Float64},d_max::Int64,iter::Int64,bf::Int64)::Tuple{Array{Float64,2},Array{Float64,2},Int64,Int64}
    terminals::Vector{Int64} = find(x->(x==1),tree[:,6])
    variables::Vector{Int64} = StatsBase.sample(1:K,Int64(round(p*K)),replace = false)
    test::Array{Tuple{Int64,Int64},2} = collect(Base.product(variables,terminals))
    gammai::Float64 = gamma[rand(1:length(gamma),1)[1]]
    fit::Dict{Float64,Vector{Float64}} = Dict{Float64,Vector{Float64}}()
    sizehint!(fit,length(test))

    terminal::Int64 = 0
    variable::Int64 = 0
    Nt::Int64 = min(20,N)
    deep::Int64 = 0
    xtest0::Vector{Float64} = Vector{Float64}(N)
    xtest::Vector{Float64} = Vector{Float64}(Nt)
    pr::Vector{Float64} = Vector{Float64}(N)
    gammascale::Float64 = 0.0
    ssr::Matrix{Float64} =  Matrix{Float64}(Nt,3+iter)
    bestres::Int64 = 1
    res0::Vector{Float64} = Vector{Float64}(5+iter)
    middlenodes::Vector{Int64} =  find(x->(isnan(x)),pmat[1,:])
    for i in 1:length(test)
        terminal = test[i][2]
        variable = test[i][1]
        deep = tree[terminal,9]
        xtest0 = x[:,variable]
        pr = pmat[:,terminal]
        xtest =  StatsBase.wsample(xtest0, pr+0.01, min(20,N),replace = false)
        gammascale = max(std(xtest0),0.1)
        for j in 1:length(xtest)
            ssr[j,:] = nvt(xtest0,y,xtest[j],gammai/gammascale,node_obs,pmat,deep+1,terminal,middlenodes)::Vector{Float64}
        end
        bestres = find(x->(x==minimum(ssr[:,1])),ssr[:,1])[1]
        res = vcat(xtest[bestres],gammai/gammascale,ssr[bestres,:])
        fit[i] = res
    end
    
    # = 3 e o best val = #
    fits::Vector{Float64} = [fit[i][3] for i in 1:length(fit)]
    best::Int64 = find(x->(x==minimum(fits)),fits)[1]
    node::Vector{Float64} = fit[best]

    if bf==5
        iter = d_max +1
        return tree, pmat, bf, iter
    end
    if node[3]==Inf
        bf = bf + 1
        return tree, pmat, bf, iter
    end

    bterminal::Int64 = test[best][2]
    
    ### 1-side, 2-b, 3-c0, 4-gamma, 5-parent, 6-terminal,7-variable, 8-id, 9-deep
    nodes::Matrix{Float64}=Matrix{Float64}(2,9) 
    nodes[1,:] = [1, node[4], node[1], node[2], tree[bterminal,8], 1, test[best][1], size(tree,1)+1, tree[bterminal,9]+1]
    nodes[2,:] = [2, node[5], node[1], node[2], tree[bterminal,8], 1, test[best][1], size(tree,1)+2, tree[bterminal,9]+1]

    # adjust old tree terminal nodes and updated bs
    tree[bterminal,6] = 0
    rem::Vector{Int64} = find(x->x==0,tree[:,6])
    terminalst::Vector{Bool} = trues(size(tree,1))
    terminalst[rem] = false
    tree[terminalst,2] = node[6:end]

    tree = vcat(tree,nodes)
    
    xb::Vector{Float64} = x[:,test[best][1]]
    probs::Vector{Float64}=Vector{Float64}(N) 
    [probs[i] = 1/(1+exp(-nodes[1,4]*(xb[i]-nodes[1,3]))) for i in 1:length(xb)]
    pmataux::Matrix{Float64} = hcat(probs,1-probs)
    [pmataux[i,j] = pmataux[i,j]*pmat[i,bterminal] for i in 1:length(xb), j in 1:2]

    pmat::Matrix{Float64} = hcat(pmat,pmataux)::Matrix{Float64}
    [pmat[i,bterminal] = NaN for i in 1:N]
    iter = iter + 1
    return tree, pmat, bf, iter
end

function grow0(x::Matrix{Float64},y::Vector{Float64},p::Float64,N::Int64,K::Int64,gamma::Vector{Float64},node_obs::Int64)::Tuple{Array{Float64,2},Array{Float64,2}}
    variables::Vector{Int64} = StatsBase.sample(1:K,Int64(round(p*K)),replace = false)
    gammai::Float64 = gamma[rand(1:length(gamma),1)[1]]
    fit0::Dict{Float64,Vector{Float64}} = Dict{Float64,Vector{Float64}}()
    sizehint!(fit0,length(variables))
    Nt::Int64 = min(20,N)
    ssr::Matrix{Float64} =  Matrix{Float64}(Nt,3)
    xtest0::Vector{Float64} = Vector{Float64}(N)
    xtest::Vector{Float64} = Vector{Float64}(Nt)
    gammascale::Float64 = 1.0::Float64
    bestres::Int64 = 1::Int64
    res0::Vector{Float64} = Vector{Float64}(5)

    for i in 1:length(variables)
        xtest0 = x[:,variables[i]]
        xtest = StatsBase.sample(xtest0, Nt,replace = false)
        gammascale = max(std(xtest0),0.1)
        for j in 1:length(xtest)
            ssr[j,:] = invt(xtest0,y,xtest[j],gammai/gammascale,node_obs)
        end

        bestres= find(x->(x==minimum(ssr[:,1])),ssr[:,1])[1]
        res0 = vcat(xtest[bestres],gammai/gammascale,ssr[bestres,:])
        fit0[i] = res0
    end

    # = 3 e o best val = #
    fit0s::Vector{Float64} = [fit0[i][3] for i in 1:length(fit0)]
    best0::Int64 = find(x->(x==minimum(fit0s)),fit0s)[1]
    node0::Vector{Float64} = fit0[best0]

    ### 1-side, 2-b, 3-c0, 4-gamma, 5-parent, 6-terminal,7-variable, 8-id, 9-deep
    node0left::Vector{Float64} = [1, node0[4], node0[1], node0[2], 0, 1, variables[best0], 1, 1]
    node0right::Vector{Float64} = [2, node0[5], node0[1], node0[2], 0, 1, variables[best0], 2, 1]

    #tree::Matrix{Float64} = (hcat(node0left,node0right))'
    tree::Matrix{Float64} = Matrix{Float64}(2,9)
    tree[1,:] = node0left
    tree[2,:] = node0right

    xb::Vector{Float64} = x[:,variables[best0]]
    probs::Vector{Float64}=Vector{Float64}(N) 
    [probs[i] = 1/(1+exp(-node0left[4]*(xb[i]-node0left[3]))) for i in 1:length(xb)]
    pmat::Matrix{Float64} = hcat(probs,1-probs)
    
    return tree, pmat
end    

function nvt(x::Vector{Float64},y::Vector{Float64},c0::Float64,gamma::Float64,node_obs::Int64,pmat::Matrix{Float64},deep::Int64,terminal::Int64,middlenodes::Vector{Int64})::Vector{Float64}
    N::Int64 = length(x)
    logit::Vector{Float64} = Vector{Float64}(N)
    [logit[i] = 1/(1+exp(-gamma*(x[i]-c0))) for i in 1:N]

    terminal_logit::Vector{Float64} = pmat[:,terminal]
    bs::Matrix{Float64} = Matrix{Float64}(N,2)
    for i in 1:N
        bs[i,1] = logit[i] .* terminal_logit[i]
        bs[i,2] = (1-logit[i]) .* terminal_logit[i]
    end
    l0::Float64 = sum(sign.(bs[:,1]-0.5^deep)+1)/2
    l1::Float64 = sum(sign.(bs[:,2]-0.5^deep)+1)/2
    
    rem::Vector{Int64} = vcat(terminal,middlenodes)
    Bools::Vector{Bool} = trues(size(pmat,2))
    Bools[rem] = false
    X::Matrix{Float64} = hcat(bs,pmat[:,Bools])

    if l0 < node_obs || l1 < node_obs
        e::Vector{Float64} = fill(Inf64,size(X,2)+1)
        return e
    end

    res::Vector{Float64} = try 
        fastlm(X,y)::Vector{Float64}
    catch
        fill(Inf64,size(X,2)+1)::Vector{Float64}
    end
    
    return res
end

function invt(x::Vector{Float64},y::Vector{Float64},c0::Float64,gamma::Float64,node_obs::Int64)::Vector{Float64}
    N::Int64 = length(x)
    b0 = Vector{Float64}(N)
    [b0[i] = 1/(1+exp(-gamma*(x[i]-c0))) for i in 1:N]
    b1::Vector{Float64} = 1 - b0

    l0::Float64 =  sum(sign.(b0-0.5)+1)/2
    l1::Float64 =  sum(sign.(b1-0.5)+1)/2
    #l0 = length(find(x->(x>=0.5),b0))
    #l1 = length(find(x->(x>=0.5),b1))

    X = Matrix{Float64}(N,2)
    X[:,1] = b0
    X[:,2] = b1

    if l0 < node_obs || l1 < node_obs
        e::Vector{Float64} = fill(Inf64,size(X,2)+1)
        return e
    end

    res::Vector{Float64} = try 
        fastlm(X,y)
    catch
        fill(Inf64,size(X,2)+1)
    end
    
    return res
end

function fastlm1(x::Matrix{Float64},y::Vector{Float64})::Vector{Float64}
    bhat::Vector{Float64} = inv(x'*x)*x'*y
    error::Vector{Float64} = y-x*bhat
    ssr::Float64 = sum( [error[i]^2 for i in 1:length(y)])
    res::Vector{Float64} = vcat(ssr,bhat)
    return res
end
