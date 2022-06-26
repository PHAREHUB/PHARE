class GridParams:
    def __init__(self,dim,interpOrder):
        self.dim = dim
        self.interpOrder = interpOrder
        self.nbrCell = ()
        self.dl = ()


    def setNbrCell(self,nbrCellX,nbrCellY,nbrCellZ):
        if (self.dim == 1):
            self.nbrCell = (nbrCellX)
        elif (self.dim == 2):
            self.nbrCell = (nbrCellX, nbrCellY)
        elif (self.dim == 3):
            self.nbrCell = (nbrCellX,nbrCellY,nbrCellZ)

    def setDl(self,dx,dy,dz):
        if (self.dim == 1):
            self.dl = (dx)
        elif (self.dim == 2):
            self.dl = (dx,dy)
        elif (self.dim == 3):
            self.dl = (dx,dy,dz)
