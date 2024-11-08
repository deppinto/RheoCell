import numpy as np
from math import atan2

def _charge_array(Type='tetratic'):
    """
    Compute the charge array of vector field (vx,vy)

    Args:
        vx: x component of vector field
        vy: y component of vector field
        Type: nematic -- (vx,vy) is head-tail symmetrical
              polar-- (vx,vy) is not head-tail symmetical
    Returns:
        Field of the charge distribution with the same shape as vx and vy
    """
    # compute angle
    def wang(a, b):
        """Infamous chinese function"""
        '''Type = 'nematic' or 'polar' '''
        ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])
        # if (Type == 'nematic') and (ang > pi/2.):
        #     b = [-i for i in b]
        while(Type == 'tetratic' and (abs(ang) > np.pi/4)):
            bx = -b[1]
            by = b[0]
            b[0] = bx
            b[1] = by
            # print(b)
            ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])
            # print(ang/np.pi)
        m = a[0]*b[1]-a[1]*b[0]
        return -np.sign(m)*atan2(abs(m), a[0]*b[0]+a[1]*b[1])

    # This mysterious part was stolen from Amin's code.
    # (calculate the winding number)
    ax1 = [1,0]
    ax2 = [np.sqrt(3),1]
    ax3 = [1,0]
    ax4 = [1,np.sqrt(3)]
    ax5 = [1,0]
    ax6 = [1,0]
    ax7 = [1,0]
    ax8 = [1,1]

    # ax1 = [1,0]
    # ax2 = [1,np.sqrt(3)]
    # ax3 = [1,0]
    # ax4 = [np.sqrt(3),1]
    # ax5 = [1,0]
    # ax6 = [1,0]
    # ax7 = [1,0]
    # ax8 = [1,1]

    w = wang(ax1, ax5)
    w += wang(ax5, ax3)
    w += wang(ax3, ax6)
    w += wang(ax6, ax2)
    w += wang(ax2, ax8)
    w += wang(ax8, ax4)
    w += wang(ax4, ax7)
    w += wang(ax7, ax1)
    w /= 2.*np.pi

    return w

print(_charge_array('tetratic'))
