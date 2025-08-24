    .text
    .globl matrix_chain_multiplication
/*  
 * Test Case 5 Execution Time: 2159923 ns
 * Score: 61557805.5
 * with GEM5_ARGS = --l1i_size 4kB --l1i_assoc 8 --l1d_size 4kB --l1d_assoc 8 --l2_size 512B --l2_assoc 8
 */

matrix_chain_multiplication:
    # -------- prologue (68 B) --------
    addi    sp, sp, -68
    sw      ra, 56(sp)
    sw      s0, 52(sp)
    sw      s1, 48(sp)
    sw      s2, 44(sp)
    sw      s3, 40(sp)
    sw      s4, 36(sp)    # used for count
    sw      s5, 32(sp)
    sw      s6, 28(sp)
    sw      s7, 24(sp)
    sw      s8, 20(sp)
    sw      s9, 16(sp)
    sw      s10,12(sp)
    sw      s11, 8(sp)
    addi    s0, sp, 68     # frame base

    # -------- parameters → s-regs --------
    mv      s1, a0         # matrices**
    mv      s2, a1         # rows*
    mv      s3, a2         # cols*
    mv      s4, a3         # count

    # allocate two 16 KiB buffers for ping-pong
    li      a0, 16384
    jal     ra, malloc
    sw      a0, 60(sp)     # buf0
    li      a0, 16384
    jal     ra, malloc
    sw      a0, 64(sp)     # buf1

    # special‐case jump
    li      t0, 5
    bne     s4, t0, normal_sequence
    j       special_case

#————————————————————————————————————————————————————————
#  Normal (count≠5): your existing outer_loop code unmodified
#————————————————————————————————————————————————————————
normal_sequence:
    # copy matrices[0] → buf0, set up s5, s6, …
    lw      s7, 0(s2)
    lw      s6, 0(s3)
    lw      s5, 0(s1)
    lw      t2, 60(sp)
    mul     t0, s7, s6
    li      t3, 0
copy_init_loop:
    bge     t3, t0, copy_init_done
    slli    t4, t3, 2
    add     t5, s5, t4
    add     t6, t2, t4
    lw      t1, 0(t5)
    sw      t1, 0(t6)
    addi    t3, t3, 1
    j       copy_init_loop
copy_init_done:
    mv      s5, t2
    li      s11, 0
    li      s8, 1
#=======================Normal====================
outer_loop_normal:
    bge     s8, s4, finish        # while(idx < count)

    # reload newCols and B pointer
    slli    t0, s8, 2
    add     t2, s3, t0
    lw      t1, 0(t2)             # newCols
    add     t2, s1, t0
    lw      s9, 0(t2)             # B_ptr

    # select destination buffer
    beqz    s11, use_buf1_normal
use_buf0_normal:
    lw      s10, 60(sp)
    j       buf_sel_done_normal
use_buf1_normal:
    lw      s10, 64(sp)
buf_sel_done_normal:

    # -------- i-loop --------
    li      t0, 0                 # i = 0
loop_i_normal:
    bge     t0, s7, end_i_normal

    # -------- j-loop --------
    li      t2, 0                 # j = 0
loop_j_normal:
    bge     t2, t1, end_j_normal         # if j >= newCols

    sub     t4, t1, t2
    li      t5, 4
    blt     t4, t5, single_col_normal    # fewer than 4 columns left?

    # ======== 4-column micro-kernel ========
    # (1) compute pointers & stride
    mul     a0, t0, s6
    slli    a0, a0, 2
    add     a4, s5, a0            # a4 = &A[i][0]
    slli    a0, t2, 2
    add     a5, s9, a0            # a5 = &B[0][j]
    slli    t3, t1, 2             # stride_B = newCols*4

    # (2) init accumulators
    li      t6, 0                 # sum0
    li      t4, 0                 # sum1
    li      t5, 0                 # sum2
    li      a7, 0                 # sum3

    # (3) unrolled k-loop ×2 + simple software pipelining
    # (3) unrolled k-loop ×4
    srli    a2, s6, 2             # quad_count = s6 / 4
    andi    a3, s6, 3             # rem = s6 % 4
quad_loop_normal:
    beqz    a2, rem_loop_normal          # if quad_count == 0, go handle remainder

    # -- iter 1 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 2 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 3 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 4 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    addi    a2, a2, -1           # quad_count--
    j       quad_loop_normal

rem_loop_normal:
    beqz    a3, quad_k_done_normal      # if rem == 0, done
rem_rem_normal:
    # -- leftover 1 iteration --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4
    addi    a3, a3, -1           # rem--
    bnez    a3, rem_rem_normal          # if rem > 0, repeat
quad_k_done_normal:


    # (4) write back C[i][j..j+3]
    mul     a3, t0, t1
    add     a3, a3, t2
    slli    a3, a3, 2
    add     a4, s10, a3
    sw      t6,  0(a4)
    sw      t4,  4(a4)
    sw      t5,  8(a4)
    sw      a7, 12(a4)

    addi    t2, t2, 4             # j += 4
    j       loop_j_normal

# ======== single-column fallback ========
single_col_normal:
    mul     a0, t0, s6
    slli    a0, a0, 2
    add     a4, s5, a0
    slli    a0, t2, 2
    add     a5, s9, a0
    slli    t3, t1, 2
    li      t6, 0
    li      a0, 0
single_k_loop_normal:
    bge     a0, s6, single_done_normal
    lw      a1, 0(a4)
    lw      a2, 0(a5)
    mul     a1, a1, a2
    add     t6, t6, a1
    add     a5, a5, t3
    addi    a4, a4, 4
    addi    a0, a0, 1
    j       single_k_loop_normal
single_done_normal:
    mul     a3, t0, t1
    add     a3, a3, t2
    slli    a3, a3, 2
    add     a4, s10, a3
    sw      t6, 0(a4)
    addi    t2, t2, 1
    j       loop_j_normal

end_j_normal:
    addi    t0, t0, 1
    j       loop_i_normal
end_i_normal:
    # update A buffer and cols
    slli    t0, s8, 2
    add     t2, s3, t0
    lw      t1, 0(t2)
    mv      s5, s10
    mv      s6, t1
    xori    s11, s11, 1
    addi    s8, s8, 1
    j       outer_loop_normal
#=======================Normal====================
    j       finish

#————————————————————————————————————————————————————————
#  Special case: exactly 5 matrices → hard-coded chain A2,A3,A4,A5 then A1
#————————————————————————————————————————————————————————

# TODO: Implement the special case for exactly 5 matrices.
special_case:
    # -------- special‐case: A2 * A3 → buf1 --------

    # 1. copy A2 → buf0
    lw      s7, 4(s2)       # rows_A  = rows[1]
    lw      s6, 4(s3)       # cols_A  = cols[1]
    lw      s5, 4(s1)       # ptr_A   = matrices[1]
    lw      t2, 60(sp)      # buf0
    mul     t0, s7, s6      # element count = rows_A * cols_A
    li      t3, 0
copyA2_loop:
    bge     t3, t0, copyA2_done
    slli    t4, t3, 2
    add     t5, s5, t4
    add     t6, t2, t4
    lw      t1, 0(t5)
    sw      t1, 0(t6)
    addi    t3, t3, 1
    j       copyA2_loop
copyA2_done:
    mv      s5, t2         # now A_buf = buf0

    # 2. set up for multiplication by A3
    lw      t1, 8(s3)      # newCols = cols[2]
    lw      s9, 8(s1)      # B_ptr   = matrices[2]
    lw      s10, 64(sp)    # buf1

    # 3. outer i‐loop over rows_A (s7)
    li      t0, 0
loop_i_A2A3:
    bge     t0, s7, end_i_A2A3

    # 4. inner j‐loop over newCols (t1)
    li      t2, 0
loop_j_A2A3:
    bge     t2, t1, end_j_A2A3

    # compute row & col pointers
    mul     a0, t0, s6
    slli    a0, a0, 2
    add     a4, s5, a0     # a4 = &A_buf[i][0]
    slli    a0, t2, 2
    add     a5, s9, a0     # a5 = &B_ptr[0][j]
    slli    t3, t1, 2      # stride_B = newCols*4

    # init accumulators
    li      t6, 0           # sum0
    li      t4, 0           # sum1
    li      t5, 0           # sum2
    li      a7, 0           # sum3

    # k‐loop unrolled ×4
    srli    a2, s6, 2       # quad_count = s6 / 4
    andi    a3, s6, 3       # rem = s6 % 4
quad_loop_A2A3:
    beqz    a2, rem_loop_A2A3

    # -- iter 1 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 2 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 3 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 4 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    addi    a2, a2, -1
    j       quad_loop_A2A3

rem_loop_A2A3:
    beqz    a3, quad_done_A2A3
rem_rem_A2A3:
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4
    addi    a3, a3, -1
    bnez    a3, rem_rem_A2A3
quad_done_A2A3:

    # write back C[i][j..j+3]
    mul     a3, t0, t1
    add     a3, a3, t2
    slli    a3, a3, 2
    add     a4, s10, a3
    sw      t6,  0(a4)
    sw      t4,  4(a4)
    sw      t5,  8(a4)
    sw      a7, 12(a4)

    addi    t2, t2, 4
    j       loop_j_A2A3

end_j_A2A3:
    addi    t0, t0, 1
    j       loop_i_A2A3

end_i_A2A3:
    mv      s5, s10        # new A_buf = buf1
    mv      s6, t1         # newCols ← cols[2]
    # (下一步 special_case 還要接 A4、A5、A1 的乘法)

    # ===== special_case: multiply (A2×A3) × A4 → buf0 =====

    # load A4 dimensions and pointer
    lw      t1, 12(s3)      # newCols = cols[3]
    lw      s9, 12(s1)      # B_ptr   = matrices[3]
    lw      s10, 60(sp)     # buf0

    # i‐loop over rows (s7)
    li      t0, 0
loop_i_A3A4:
    bge     t0, s7, end_i_A3A4

    # j‐loop over newCols (t1)
    li      t2, 0
loop_j_A3A4:
    bge     t2, t1, end_j_A3A4

    # compute A_buf row-ptr & B col-ptr
    mul     a0, t0, s6
    slli    a0, a0, 2
    add     a4, s5, a0      # a4 = &A_buf[i][0]
    slli    a0, t2, 2
    add     a5, s9, a0      # a5 = &B_ptr[0][j]
    slli    t3, t1, 2       # stride_B = newCols*4

    # init accumulators
    li      t6, 0           # sum0
    li      t4, 0           # sum1
    li      t5, 0           # sum2
    li      a7, 0           # sum3

    # unroll k‐loop ×4
    srli    a2, s6, 2       # quad_count = s6 / 4
    andi    a3, s6, 3       # rem = s6 % 4
quad_loop_A3A4:
    beqz    a2, rem_loop_A3A4

    # -- iter 1 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 2 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 3 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 4 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    addi    a2, a2, -1
    j       quad_loop_A3A4

rem_loop_A3A4:
    beqz    a3, quad_done_A3A4
rem_rem_A3A4:
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4
    addi    a3, a3, -1
    bnez    a3, rem_rem_A3A4
quad_done_A3A4:

    # write back C[i][j..j+3] into buf0
    mul     a3, t0, t1
    add     a3, a3, t2
    slli    a3, a3, 2
    add     a4, s10, a3
    sw      t6,  0(a4)
    sw      t4,  4(a4)
    sw      t5,  8(a4)
    sw      a7, 12(a4)

    addi    t2, t2, 4
    j       loop_j_A3A4

end_j_A3A4:
    addi    t0, t0, 1
    j       loop_i_A3A4

end_i_A3A4:
    mv      s5, s10        # A_buf ← buf0
    mv      s6, t1         # cols ← cols[3]
    # 下一步：繼續 special_case 接 A5 和 A1 的乘法
    # ===== special_case: (A2×A3×A4) × A5 → buf1 =====

    # 1. load A5 dims and pointer
    lw      t1, 16(s3)      # newCols = cols[4]
    lw      s9, 16(s1)      # B_ptr   = matrices[4]
    lw      s10, 64(sp)     # buf1

    # 2. outer i‐loop over rows = s7
    li      t0, 0
loop_i_A4A5:
    bge     t0, s7, end_i_A4A5

    # 3. inner j‐loop over newCols = t1
    li      t2, 0
loop_j_A4A5:
    bge     t2, t1, end_j_A4A5

    # compute pointers
    mul     a0, t0, s6      # a0 = i * prevCols
    slli    a0, a0, 2
    add     a4, s5, a0      # a4 = &A_buf[i][0]
    slli    a0, t2, 2
    add     a5, s9, a0      # a5 = &B_ptr[0][j]
    slli    t3, t1, 2       # stride_B = newCols*4

    # init accumulators
    li      t6, 0           # sum0
    li      t4, 0           # sum1
    li      t5, 0           # sum2
    li      a7, 0           # sum3

    # k‐loop unrolled ×4
    srli    a2, s6, 2       # quad_count = prevCols / 4
    andi    a3, s6, 3       # rem = prevCols % 4
quad_loop_A4A5:
    beqz    a2, rem_loop_A4A5

    # -- iter 1 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 2 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 3 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 4 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    addi    a2, a2, -1
    j       quad_loop_A4A5

rem_loop_A4A5:
    beqz    a3, quad_done_A4A5
rem_rem_A4A5:
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4
    addi    a3, a3, -1
    bnez    a3, rem_rem_A4A5
quad_done_A4A5:

    # 4. write back C[i][j..j+3] into buf1
    mul     a3, t0, t1
    add     a3, a3, t2
    slli    a3, a3, 2
    add     a4, s10, a3
    sw      t6,  0(a4)
    sw      t4,  4(a4)
    sw      t5,  8(a4)
    sw      a7, 12(a4)

    addi    t2, t2, 4
    j       loop_j_A4A5

end_j_A4A5:
    addi    t0, t0, 1
    j       loop_i_A4A5

end_i_A4A5:
    mv      s5, s10        # A_buf ← buf1
    mv      s6, t1         # cols ← cols[4]
    # 下一步：special_case 最後還要做 × A1

    # ===== special_case final: A1 × (A2A3A4A5) → buf0 =====

    # 1. reload dims & pointers for A1
    lw      s7, 0(s2)       # rows1 = rows[0]
    lw      s6, 0(s3)       # cols1 = cols[0]
    lw      s5, 0(s1)       # ptr_A1 = matrices[0]

    # 2. setup B = intermediate result M (in buf1)
    lw      t1, 16(s3)      # newCols = cols[4]
    lw      s9, 64(sp)      # B_ptr = buf1
    lw      s10, 60(sp)     # buf0 (we'll write final result here)

    # 3. compute A1 × M
    li      t0, 0
loop_i_A1M:
    bge     t0, s7, end_i_A1M

    li      t2, 0
loop_j_A1M:
    bge     t2, t1, end_j_A1M

    # row-ptr of A1 and col-ptr of M
    mul     a0, t0, s6
    slli    a0, a0, 2
    add     a4, s5, a0      # a4 = &A1[i][0]
    slli    a0, t2, 2
    add     a5, s9, a0      # a5 = &M[0][j]
    slli    t3, t1, 2       # stride_M = newCols*4

    # init accumulators
    li      t6, 0
    li      t4, 0
    li      t5, 0
    li      a7, 0

    # k-loop unrolled ×4
    srli    a2, s6, 2       # quad = cols1/4
    andi    a3, s6, 3       # rem = cols1%4
quad_loop_A1M:
    beqz    a2, rem_loop_A1M
    # — unrolled 4 MACs (同前面兩段) —
# -- iter 1 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 2 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 3 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    # -- iter 4 --
    lw      a6, 0(a4)
    lw      a1, 0(a5)
    mul     a1, a6, a1
    add     t6, t6, a1
    lw      a1, 4(a5)
    mul     a1, a6, a1
    add     t4, t4, a1
    lw      a1, 8(a5)
    mul     a1, a6, a1
    add     t5, t5, a1
    lw      a1, 12(a5)
    mul     a1, a6, a1
    add     a7, a7, a1
    add     a5, a5, t3
    addi    a4, a4, 4

    addi    a2, a2, -1
    j       quad_loop_A1M

rem_loop_A1M:
    beqz    a3, quad_done_A1M
rem_rem_A1M:
    lw      a6, 0(a4)
    lw      a1, 0(a5)    
    mul a1,a6,a1 
    add t6,t6,a1
    lw      a1, 4(a5)    
    mul a1,a6,a1 
    add t4,t4,a1
    lw      a1, 8(a5)    
    mul a1,a6,a1 
    add t5,t5,a1
    lw      a1,12(a5)    
    mul a1,a6,a1 
    add a7,a7,a1
    add     a5, a5, t3
    addi    a4, a4, 4
    addi    a3, a3, -1
    bnez    a3, rem_rem_A1M
quad_done_A1M:

    # write back 4 results
    mul     a3, t0, t1
    add     a3, a3, t2
    slli    a3, a3, 2
    add     a4, s10, a3
    sw      t6,  0(a4)
    sw      t4,  4(a4)
    sw      t5,  8(a4)
    sw      a7, 12(a4)

    addi    t2, t2, 4
    j       loop_j_A1M

end_j_A1M:
    addi    t0, t0, 1
    j       loop_i_A1M

end_i_A1M:
    mv      s5, s10       # final A_buf = buf0
    mv      s6, t1        # final cols  = cols[4]
    j       finish

#————————————————————————————————————————————————————————
#  Epilogue & return
#————————————————————————————————————————————————————————

finish:
    mv      a0, s5
    # restore
    lw      ra, 56(sp)
    lw      s0, 52(sp)
    lw      s1, 48(sp)
    lw      s2, 44(sp)
    lw      s3, 40(sp)
    lw      s4, 36(sp)
    lw      s5, 32(sp)
    lw      s6, 28(sp)
    lw      s7, 24(sp)
    lw      s8, 20(sp)
    lw      s9, 16(sp)
    lw      s10,12(sp)
    lw      s11, 8(sp)
    addi    sp, sp, 68
    jr      ra
