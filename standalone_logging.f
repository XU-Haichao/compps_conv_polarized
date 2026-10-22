      subroutine xwrite(message,level)
      implicit none
      character*(*) message
      integer level
c Standalone ctypes builds do not link against XSPEC's logging service.
      return
      end

      subroutine fxwrite(message,level)
      implicit none
      character*(*) message
      integer level
      return
      end
