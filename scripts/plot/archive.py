# This file is part of CELADRO, Copyright (C) 2016-17, Romain Mueller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import archive_base


class archive(archive_base.archive):
    """Simply reshape 2d fields after importing"""

    def __init__(self, path):
        archive_base.archive.__init__(self, path)
        self.parameters['walls'] = np.reshape(self.parameters['walls'],
                                              self.Size)
        self.__dict__.update(self.parameters)

    def unfold_patch(self,frame,name):
        def getElementX(site, distX, sidex, sidey):
            x = int( (site-(int(site/sidex)*sidex))+distX )
            if x<0:
             x+=sidex
            if x>=sidex:
                x-=sidex
            return x

        def getElementY(site, distY, sidex, sidey):
            y = int( int(site/sidex)+distY )
            if y<0:
                y+=sidey
            if y>=sidey:
                y-=sidey
            return y

        def getElement(site, LsubX, LsubY, offsetX, offsetY, sub_corner_bottom_left, boxXsize, boxYsize):
            return getElementX(sub_corner_bottom_left, ((site - int(site/LsubX) * LsubX)+offsetX)%LsubX, boxXsize, boxYsize ) + getElementY(sub_corner_bottom_left, ((site/LsubX)+offsetY)%LsubY , boxXsize, boxYsize) * boxXsize

        def getX(site, LsubX, LsubY, offsetX, offsetY, sub_corner_bottom_left, boxXsize, boxYsize):
            return getElementX(sub_corner_bottom_left, ((site - int(site/LsubX) * LsubX)+offsetX)%LsubX , boxXsize, boxYsize)

        def getY(site, LsubX, LsubY, offsetX, offsetY, sub_corner_bottom_left, boxXsize, boxYsize):
            return getElementY(sub_corner_bottom_left, ((site/LsubX)+offsetY)%LsubY , boxXsize, boxYsize)
        

        tmp_patch = getattr(frame,name)
        lx, ly = self.Size
        #px, py = self.patch_size
        rtn = []
        for i in range(len(tmp_patch)):
            px = self.patch_size[i][0]
            py = self.patch_size[i][1]
            p=[[0. for q in range(lx)] for k in range(ly)]
            patch=[]
            for j in range(len(tmp_patch[i])):
                patch.append(float(tmp_patch[i][j]))

            cornerSite=int(frame.patch_min[i][0])+int(frame.patch_min[i][1])*lx
            for j in range(py):
                for k in range(px):
                    row_site=k+j*px
                    x=getX(row_site, px, py, int(frame.offset[i][0]), int(frame.offset[i][1]), cornerSite, lx, ly)
                    y=getY(row_site, px, py, int(frame.offset[i][0]), int(frame.offset[i][1]), cornerSite, lx, ly)
                    p[y][x]=patch[row_site]


            '''
            p = np.reshape(patch, (py, px))
            # compensate for offset
            p = np.roll(p, frame.offset[i][0], axis=0)
            p = np.roll(p, frame.offset[i][1], axis=1)
            # extend to full size
            p = np.concatenate((p, np.zeros((px, ly-py))), axis=0)
            p = np.concatenate((p, np.zeros((lx-px, ly))), axis=1)
            # put in right postition
            p = np.roll(p, frame.patch_min[i][0], axis=0)
            p = np.roll(p, frame.patch_min[i][1], axis=1)
            # save
            '''
            rtn.append(p)

        '''
        cc=0
        for i in range(len(rtn[1])):
            if rtn[1][i]>0:
                print(cc, rtn[1][i])
            cc+=1
        '''

        return rtn
        

    def read_frame(self, frame):
        frame = archive_base.archive.read_frame(self, frame)

        # array sizes
        lx, ly = self.Size
        #px, py = self.patch_size

        '''
        phi = []
        for i in range(len(frame.phi)):
            p = np.reshape(frame.phi[i], (px, py))
            # compensate for offset
            p = np.roll(p, frame.offset[i][0], axis=0)
            p = np.roll(p, frame.offset[i][1], axis=1)
            # extend to full size
            p = np.concatenate((p, np.zeros((px, ly-py))), axis=1)
            p = np.concatenate((p, np.zeros((lx-px, ly))), axis=0)
            # put in right postition
            p = np.roll(p, frame.patch_min[i][0], axis=0)
            p = np.roll(p, frame.patch_min[i][1], axis=1)
            # save
            phi.append(p)
        '''

        frame.phi = self.unfold_patch(frame,'phi')
        for var_name in ['fp_x','fp_y','fpol_x','fpol_y','fnem_x','fnem_y','fshape_x','fshape_y']:
            if hasattr(frame,var_name):
                setattr(frame,var_name,self.unfold_patch(frame,var_name))


        #reshape the velocity field and force density
        if hasattr(frame,'fdipole_field_x'):
            frame.fdipole_field_x = np.reshape(frame.fdipole_field_x,(lx,ly))
            frame.fdipole_field_y = np.reshape(frame.fdipole_field_y,(lx,ly))
        if hasattr(frame,'fpol_field_x'):
            frame.fpol_field_x = np.reshape(frame.fpol_field_x,(lx,ly))
            frame.fpol_field_y = np.reshape(frame.fpol_field_y,(lx,ly))
        if hasattr(frame,'fp_field_x'):
            frame.fp_field_x = np.reshape(frame.fp_field_x,(lx,ly))
            frame.fp_field_y = np.reshape(frame.fp_field_y,(lx,ly))
        if hasattr(frame,'velocity_field_x'):
            frame.velocity_field_x = np.reshape(frame.velocity_field_x,(lx,ly))
            frame.velocity_field_y = np.reshape(frame.velocity_field_y,(lx,ly))
        return frame


def loadarchive(path):
    return archive(path)
