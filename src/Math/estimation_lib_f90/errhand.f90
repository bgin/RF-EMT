      subroutine errhand (message);
      character(*),intent(in) :: message;
      write (6,'(a)') trim(message);
      write (*,'(a)') trim(message);
      stop "error";
      end subroutine errhand;
