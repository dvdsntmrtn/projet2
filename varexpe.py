class Expe(object):

    def __init__(self, project, model, name, member=1, ybeg=1850, yend=1850, expe_control=None, color='k', marker='.', linestyle='-'):
        self.project = project # CLIMAF project name (pre-existing or user)
        self.model = model
        self.name = name
        self.member = member
        self.ybeg = ybeg
        self.yend = yend
        self.expe_control = expe_control
        self.color = color
        self.marker = marker
        self.linestyle = linestyle
        
    def __str__(self):
        print self.name, self.ybeg

    def period(self):
        return str(self.ybeg)+'-'+str(self.yend)

    def expid(self):
        return str(self.model)+'_'+str(self.name)+'_r'+str(self.member)


class Variable(object):

    def __init__(self, name=None, table=None, grid='gr'):
        self.name = name
        self.table = table
        self.grid = grid

    def varid(self):
        return str(self.name)+'_'+str(self.table)

    def __str__(self):
        print self.name, self.table
