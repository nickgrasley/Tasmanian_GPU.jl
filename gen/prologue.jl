mutable struct TasmanianSG
    pGrid      :: Ptr{Nothing}
    dimensions :: Int
    outputs    :: Int
    depth      :: Int

    function TasmanianSG(dims::Int,nout::Int,depth::Int)
	this = new()
	output_ptr = ccall(
	    (:tsgConstructTasmanianSparseGrid,TASlib), # name of C function and library
	    Ptr{TasmanianSG},                          # output type
	    ()                                         # tuple of input types
	)
	if output_ptr == C_NULL # Could not allocate memory
	    throw(OutOfMemoryError())
	else
	    this.pGrid = output_ptr
	end
        this.dimensions = dims
        this.outputs    = nout
        this.depth      = depth
	return this
    end
end



