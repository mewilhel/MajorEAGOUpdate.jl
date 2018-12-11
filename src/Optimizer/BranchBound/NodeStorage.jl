```
    DefaultStorage

Stores the two nodes to the stack.
```
function DefaultStorage(x::Optimizer,y1::NodeBB,y2::NodeBB)
    x.MaximumNodeID += 1; x.CurrentNodeCount += 1; x.Stack[x.MaximumNodeID] = y1
    x.MaximumNodeID += 1; x.CurrentNodeCount += 1; x.Stack[x.MaximumNodeID] = y2
end

```
    DefaultStorage

Stores the one nodes to the stack.
```
function SingleStorage!(x::Optimizer,y::NodeBB)
    x.MaximumNodeID += 1; x.CurrentNodeCount += 1; x.Stack[x.MaximumNodeID] = y
end
