# states used in depth first search
@enum VISIT_STATUS VISITED UNVISITED INPROGRESS

# performs a depth-first search label nodes with distance from expression
function dfs_variable(x::Tape)
    status = Dict{Int,VISIT_STATUS}()
    parent = Dict{Int,Int}(1 => -1)
    for i in 1:length(x.nd)
        status[i] = UNVISITED
    end
    dfs_variable_kernel!(x,x.nd[1],1,status,parent)
    return parent
end

# kernel for depth first search
function dfs_variable_kernel!(x::Tape, n::Tracer.NodeData, indx, status, parent)
    first_child = children(n)[1]
    if (first_child == -1) || (first_child == -2)
        return
    else
        status[indx] = INPROGRESS
        for c in children(n)
            if status[c] == UNVISITED
                parent[c] = indx
                child_node = x.nd[c]
                dfs_variable_kernel!(x, child_node, c, status, parent)
            end
        end
        status[indx] = VISITED
    end
end

# takes a tracer tape, rearranges it to have the terminal node at the first
# index, replacing the array of child with an one elemnt array containing
# parents (PASSING)
function child_to_parent!(x::Tape)

    # reverse array order
    idx_arr = [i for i in x.set_trace_count:-1:1]
    nd_rev_temp = reverse(x.nd)
    nd_rev = Tracer.NodeData[]
    for (i,node) in enumerate(nd_rev_temp)
        if node.children[1] == -2
            push!(nd_rev,node)
        elseif node.children[1] == -1
            push!(nd_rev,node)
        else
            child_temp = [idx_arr[node.children[j]] for j in 1:length(node.children)]
            nd_temp = NodeData(node.nodetype, node.index, child_temp)
            push!(nd_rev,nd_temp)
        end
    end
    x.nd = nd_rev

    num_valued = Dict{Int,Bool}()
    for i in 1:x.set_trace_count
        num_valued[i] = x.num_valued[x.set_trace_count-i+1]
    end
    x.num_valued = num_valued

    # rewrite node list to have only parents not children
    parent_dict = dfs_variable(x)
    nd_list = []
    for i in 1:length(x.nd)
        temp_node = Tracer.NodeData(x.nd[i].nodetype,x.nd[i].index,[parent_dict[i]])
        push!(nd_list, temp_node)
    end
    x.nd = nd_list
end

function load_evaluator()
end
