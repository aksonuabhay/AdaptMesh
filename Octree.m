classdef Octree < handle
   
    properties (SetAccess = private)
        NodeCount; % Total number of nodes created.
        NodeBoundaries; % NodeCount-by-6 [MIN MAX] coordinates of node edges.
        Point3Ds; % The coordinate of points in the decomposition
        Point3DNodes; % Indices of the node that each point belongs to.
        NodeDepths; % The number of subdivisions to reach each node.
        NodeParents = zeros(0,1); % Indices of the node that each node belongs to.
        Properties; % Name/Val pairs used for creation
    end
    
    methods
        
        function this = Octree(pts,varargin)
            validateattributes(pts,{'numeric'},...
                {'real','finite','nonnan','ncols', 3},...
                mfilename,'PTS')
            
            ptCount = size(pts,1);
            this.NodeBoundaries = [min(pts,[],1) max(pts,[],1)];
            this.Point3Ds = pts;
            this.Point3DNodes = ones(ptCount,1);
            this.NodeDepths = 0;
            this.NodeParents(1) = 0;
            this.NodeCount = 1;
            
            parser = inputParser;
            parser.addParamValue('nodeCapacity',ceil(ptCount)/10);
            parser.addParamValue('maxDepth',inf);
            parser.addParamValue('maxSize',inf);
            parser.addParamValue('minSize',1000 * eps);
            parser.addParamValue('style','equal');
            parser.parse(varargin{:});
            this.Properties = parser.Results;
            
            % If empty or trivial nodes return
            if ptCount<2, return; end

            this.preAllocSpace;
            this.subDivide(1);
            this.deallocSpace;
        end
        
        function deallocSpace(this)
            this.NodeBoundaries(this.NodeCount+1:end,:) = [];
            this.NodeDepths(this.NodeCount+1:end) = [];
            this.NodeParents(this.NodeCount+1:end) = [];
        end
        
        function preAllocSpace(this)
            ptCount = size(this.Point3Ds,1);
            nodesCount = ptCount;
            if isfinite(this.Properties.nodeCapacity)
                nodesCount = ceil(2*ptCount/this.Properties.nodeCapacity);
            end
            this.NodeDepths(nodesCount) = 0;
            this.NodeParents(nodesCount) = 0;
            this.NodeBoundaries(nodesCount,1) = 0;
        end
              
        function subDivide(this, begNodes)
           
            for i = 1:length(begNodes)
                nodeNo = begNodes(i);
                
                if this.NodeDepths(nodeNo)+1 >= this.Properties.maxDepth
                    continue;
                end
                              
                thisLimits = this.NodeBoundaries(nodeNo,:);
                nodeEdgeSize = diff(thisLimits([1:3;4:6]));
                maxEdgeSize = max(nodeEdgeSize);
                minEdgeSize = min(nodeEdgeSize);
                
                if minEdgeSize < this.Properties.minSize
                    continue;
                end
                
                oldCount = this.NodeCount;
                if nnz(this.Point3DNodes==nodeNo) > this.Properties.nodeCapacity
                    this.subDivideNode(nodeNo);
                    this.subDivide(oldCount+1:this.NodeCount);
                    continue;
                end
                
                if maxEdgeSize > this.Properties.maxSize
                    this.subDivideNode(nodeNo);
                    this.subDivide(oldCount+1:this.NodeCount);
                    continue;
                end
            end
        end
        
        function subDivideNode(this,nodeNo)
            nodePtMask = this.Point3DNodes==nodeNo;
            thisNodesPoint3Ds = this.Point3Ds(nodePtMask,:);
            
            oldMax = this.NodeBoundaries(nodeNo,4:6);
            oldMin = this.NodeBoundaries(nodeNo,1:3);
                        
            if strcmp('weighted',this.Properties.style) && any(nodePtMask)
                newDiv = mean(thisNodesPoint3Ds,1);
            else
                newDiv = mean([oldMin; oldMax], 1);
            end
            
            minMidMax = [oldMin newDiv oldMax];
            newBounds = minMidMax([...
                1 2 3 4 5 6;
                1 2 6 4 5 9;
                1 5 3 4 8 6;
                1 5 6 4 8 9;
                4 2 3 7 5 6;
                4 2 6 7 5 9;
                4 5 3 7 8 6;
                4 5 6 7 8 9]);
            
            nodeMap = cat(3,[0 0 0],[0 0 1],[0 1 0],[0 1 1],...
                [1 0 0],[1 0 1],[1 1 0],[1 1 1]);
            gtMask = bsxfun(@gt, thisNodesPoint3Ds, newDiv);
            [~,nodeAssignment] = max(all(bsxfun(@eq,gtMask,nodeMap),2),[],3);

            newNodeInds = this.NodeCount+1:this.NodeCount+8;
            this.NodeBoundaries(newNodeInds,:) = newBounds;
            this.NodeDepths(newNodeInds) = this.NodeDepths(nodeNo)+1;
            this.NodeParents(newNodeInds) = nodeNo;
            this.Point3DNodes(nodePtMask) = newNodeInds(nodeAssignment);
            this.NodeCount = this.NodeCount + 8;
        end
        
        function shrink(this)
            nodeChildren = arrayfun(@(i)find(this.NodeParents==i),1:this.NodeCount,'Un',0)';
            nodeIsLeaf = cellfun(@isempty, nodeChildren);
            for i = find(nodeIsLeaf(:))'
                nodeShrinkRec(i, true)
            end
            
            function nodeShrinkRec(nodeNo, isLeafNode)
                
                oldBoundaryMax = this.NodeBoundaries(nodeNo,4:6);
                oldBoundaryMin = this.NodeBoundaries(nodeNo,1:3);
                
                if isLeafNode
   
                    ptsMask = this.Point3DNodes==nodeNo;
                    if ~any(ptsMask)
                        proposedBoundaries = [oldBoundaryMin oldBoundaryMin];
                    else
                        pts = this.Point3Ds(ptsMask,:);
                        proposedBoundaries = [...
                            max([oldBoundaryMin; min(pts,[],1)]) ...
                            min([oldBoundaryMax; max(pts,[],1)])];
                    end
                    
                else
                    childBoundaries = this.NodeBoundaries(nodeChildren{nodeNo},:);
                    proposedBoundaries = [min(childBoundaries(:,1:3),[],1) max(childBoundaries(:,4:6),[],1)];
                end
                
                if ~isequal(proposedBoundaries, [oldBoundaryMin oldBoundaryMax])
                    
                    this.NodeBoundaries(nodeNo,:) = proposedBoundaries;
                    parentNode = this.NodeParents(nodeNo);
                    if parentNode>0
                        nodeShrinkRec(parentNode, false)
                    end
                    
                end
            end
        end
        
        function nodeNos = qry(this, newPts, qryDepth)

            if nargin<3
                qryDepth = max(this.NodeDepths);
            end
            
            ptCount = size(newPts,1);
            newPts = permute(newPts,[3 2 1]);
            nodeNos = ones(ptCount,1)*-1;
                        
            nodeChildren = arrayfun(@(i)find(this.NodeParents==i),1:this.NodeCount,'Un',0)';
            nodeIsLeaf = cellfun(@isempty, nodeChildren);
            ptQryRec(1:ptCount, this.NodeParents==0, 0)
            
            function ptQryRec(newIndsToCheck_, nodesToCheck, depth)
                boundsToCheck = this.NodeBoundaries(nodesToCheck,:);
                [ptInBounds, subNodeNo] = max(all(...
                    bsxfun(@ge, newPts(:,:,newIndsToCheck_), boundsToCheck(:,1:3)) & ...
                    bsxfun(@le, newPts(:,:,newIndsToCheck_), boundsToCheck(:,4:6))...
                    ,2),[],1);
            
                if ~all(ptInBounds)
                    nodeNos(newIndsToCheck_(~ptInBounds)) = -1;
                    newIndsToCheck_(~ptInBounds) = [];
                    subNodeNo(~ptInBounds) = [];
                end
                nodeNosToAssign = nodesToCheck(subNodeNo);
                newIndsToAssign = newIndsToCheck_;
                nodeNos(newIndsToAssign) = nodeNosToAssign;
                
                if depth>=qryDepth
                    return;
                end
                
                [unqNodeNos, ~, unqGrpNos] = unique(nodeNosToAssign);
                for i = 1:length(unqNodeNos)
                    thisPtMask = unqGrpNos==i;
                    if ~nodeIsLeaf(unqNodeNos(i))
                        ptQryRec(newIndsToCheck_(thisPtMask), nodeChildren{unqNodeNos(i)}, depth+1)
                    end
                end
                
            end
        end     
    end
end