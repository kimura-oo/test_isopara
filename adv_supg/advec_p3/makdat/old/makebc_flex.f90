program makebc_func
      implicit none
      double precision :: eps
      integer :: node, nelm
      double precision, allocatable :: xx(:), yy(:)
      integer :: ifbcu, sfbcu
      integer, allocatable :: nfbcu(:)
      double precision, allocatable :: ffbcu(:)
      double precision :: xmax, xmin, ymax, ymin
      integer :: i, j
      character(50) :: mesfile, bdcfile
    
    ! 関数を外部関数として宣言
      !double precision, external :: f

    ! ファイル名読み込み
      open(9, file='makebc.file', status='unknown')
      read(9, '(a)') mesfile
      read(9, '(a)') bdcfile
      close(9)

    ! メッシュファイルと境界条件出力ファイルを開く
      open(10, file=mesfile, status='old')
      open(50, file=bdcfile, status='replace')

    ! 節点数と要素数を読み込み
      read(10, *) node, nelm
      allocate(xx(node), yy(node))

    ! 全節点の座標を読み込み
      do i = 1, node
            read(10, *) j, xx(i), yy(i)
      end do
      close(10)

      eps = 1.0d-10
      sfbcu = 0

    ! 領域の最大・最小座標を取得
      xmax = maxval(xx)
      xmin = minval(xx)
      ymax = maxval(yy)
      ymin = minval(yy)

    ! 境界上の節点数をカウント
      do i = 1, node
            if ((dabs(xx(i) - xmin) .le. eps) .or. (dabs(xx(i) - xmax) .le. eps) .or. (dabs(yy(i) - ymin) .le. eps) .or. (dabs(yy(i) - ymax) .le. eps)) then
                  sfbcu = sfbcu + 1
            end if
      end do
    
    ! 境界条件を格納する配列を確保
      allocate(nfbcu(sfbcu), ffbcu(sfbcu))

      ifbcu = 0

    ! 境界上の節点を探索し、関数f(x,y)から値を設定
      do i = 1, node
            if ((dabs(xx(i) - xmin) .le. eps) .or. (dabs(xx(i) - xmax) .le. eps) .or. (dabs(yy(i) - ymin) .le. eps) .or. (dabs(yy(i) - ymax) .le. eps)) then
                  ifbcu = ifbcu + 1
                  nfbcu(ifbcu) = i
            ! ★変更点：関数f(x,y)を用いて境界条件の値を計算
                  ffbcu(ifbcu) = f(xx(i), yy(i))
            end if
      end do

    ! 境界条件ファイルに書き出し
      write(50, 120) sfbcu
      write(50, 121) (i, nfbcu(i), ffbcu(i), i = 1, sfbcu)
      close(50)

120 format (i9)
121 format (2i9, e15.6)

contains
    ! ★追加点：境界条件の値を計算する関数f(x,y)
    ! この関数内の計算式を自由に変更してください
    function f(x, y)
        implicit none
        double precision, intent(in) :: x, y
        double precision :: f

        ! 例： f(x, y) = x * y
        f = x**3 - 3 * x * y**2 + 2 * x**2 + y
    end function f

end program makebc_func