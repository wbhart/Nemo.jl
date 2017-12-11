function test_fq_nmod_embed()

    # Prelude : Creation of the finite fields and the nodes

    p = 11

    k2, x2 = FiniteField(p, 2, "x2")
    k3, x3 = FiniteField(p, 3, "x3")
    k4, x4 = FiniteField(p, 4, "x4")
    k5, x5 = FiniteField(p, 5, "x5")
    k6, x6 = FiniteField(p, 6, "x6")
    k8, x8 = FiniteField(p, 8, "x8")
    k9, x9 = FiniteField(p, 9, "x9")
    k10, x10 = FiniteField(p, 10, "x10")
    k12, x12 = FiniteField(p, 12, "x12")
    k16, x16 = FiniteField(p, 16, "x16")
    k18, x18 = FiniteField(p, 18, "x18")
    k24, x24 = FiniteField(p, 24, "x24")
    k35, x35 = FiniteField(p, 35, "x35")
    k50, x50 = FiniteField(p, 50, "x50")
    k70, x70 = FiniteField(p, 70, "x70")
    k210, x210 = FiniteField(p, 210, "x210")

    #=
    # A first test where the user ask for the embedding in the "right"
    # order

    f2_8 = embed(k2, k8)
    f2_4 = embed(k2, k4)
    f4_8 = embed(k4, k8)

    @test f2_8(x2) == f4_8(f2_4(x2))

    f2_16 = embed(k2, k16)
    f4_16 = embed(k4, k16)

    @test f2_16(x2) == f4_16(f2_4(x2))

    f8_16 = embed(k8, k16)

    @test f2_16(x2) == f8_16(f2_8(x2))
    @test f4_16(x4) == f8_16(f4_8(x4))

    f2_6 = embed(k2, k6)
    f3_6 = embed(k3, k6)
    f2_12 = embed(k2, k12)
    f3_12 = embed(k3, k12)
    f6_12 = embed(k6, k12)

    @test f2_12(x2) == f6_12(f2_6(x2))
    @test f3_12(x3) == f6_12(f3_6(x3))

    f4_12 = embed(k4, k12)

    @test f2_12(x2) == f4_12(f2_4(x2))

    f2_24 = embed(k2, k24)
    f4_24 = embed(k4, k24)
    f8_24 = embed(k8, k24)
    f3_24 = embed(k3, k24)
    f12_24 = embed(k12, k24)
    f6_12 = embed(k6, k12)
    f2_6 = embed(k2, k6)
    f2_4 = embed(k2, k4)
    f2_12 = embed(k2, k12)
    f2_24 = embed(k2, k24)
    f6_24 = embed(k6, k24)

    @test f2_24(x2) == f8_24(f2_8(x2))
    @test f2_24(x2) == f12_24(f2_12(x2))
    @test f2_24(x2) == f4_24(f2_4(x2))
    @test f4_24(x4) == f8_24(f4_8(x4))
    @test f3_24(x3) == f6_24(f3_6(x3))
    @test f4_24(x4) == f12_24(f4_12(x4))
    @test f3_24(x3) == f12_24(f3_12(x3))
    @test f6_24(x6) == f12_24(f6_12(x6))
    @test f2_24(x2) == f4_24(f2_4(x2))
    @test f2_24(x2) == f6_24(f2_6(x2))

    f2_18 = embed(k2, k18)
    f3_18 = embed(k3, k18)
    f6_18 = embed(k6, k18)
    f3_9 = embed(k3, k9)
    f9_18 = embed(k9, k18)

    @test f2_18(x2) == f6_18(f2_6(x2))
    @test f3_18(x3) == f9_18(f3_9(x3))
    @test f3_18(x3) == f6_18(f3_6(x3))

    =#

    #= "Simple" square

    f4_12 = embed(k4, k12)
    f6_12 = embed(k6, k12)
    f2_4 = embed(k2, k4)
    f2_6 = embed(k2, k6)

    @test f6_12(f2_6(x2)) == f4_12(f2_4(x2))

    =#
    
    #= A simple transitive test

    f2_4 = embed(k2, k4)
    f4_8 = embed(k4, k8)
    f8_16 = embed(k8, k16)
    f2_8 = embed(k2, k8)
    f2_16 = embed(k2, k16)
    f4_16 = embed(k4, k16)

    @test f4_8(f2_4(x2)) == f2_8(x2)
    @test f4_16(f2_4(x2)) == f2_16(x2)
    @test f8_16(f2_8(x2)) == f2_16(x2)
    @test f8_16(f4_8(x4)) == f4_16(x4)

    =#

    #= Some other triangles

    f4_8 = embed(k4, k8)
    f2_8 = embed(k2, k8)
    f2_4 = embed(k2, k4)

    f3_12 = embed(k3, k12)
    f6_12 = embed(k6, k12)
    f3_6 = embed(k3, k6)

    @test f4_8(f2_4(x2)) == f2_8(x2)
    @test f6_12(f3_6(x3)) == f3_12(x3)

    =#

    # A pretty random lattice with small numbers

    f4_12 = embed(k4, k12)
    f6_24 = embed(k6, k24)
    f8_16 = embed(k8, k16)
    f2_16 = embed(k2, k16)
    f3_6 = embed(k3, k6)
    f9_18 = embed(k9, k18)
    f6_12 = embed(k6, k12)
    f3_18 = embed(k3, k18)
    f4_16 = embed(k4, k16)
    f6_18 = embed(k6, k18)
    f2_8 = embed(k2, k8)
    f12_24 = embed(k12, k24)
    f3_9 = embed(k3, k9)
    f2_6 = embed(k2, k6)
    f8_24 = embed(k8, k24)
    f3_24 = embed(k3, k24)
    f2_12 = embed(k2, k12)
    f3_12 = embed(k3, k12)
    f2_24 = embed(k2, k24)
    f2_4 = embed(k2, k4)
    f4_24 = embed(k4, k24)
    f4_8 = embed(k4, k8)

    @test f6_18(f3_6(x3)) == f3_18(x3)
    @test f8_16(f2_8(x2)) == f2_16(x2)
    @test f12_24(f6_12(x6)) == f6_24(x6)
    @test f9_18(f3_9(x3)) == f3_18(x3)
    @test f6_24(f3_6(x3)) == f3_24(x3)
    @test f6_12(f2_6(x2)) == f2_12(x2)
    @test f6_12(f3_6(x3)) == f3_12(x3)
    @test f12_24(f3_12(x3)) == f3_24(x3)
    @test f6_24(f2_6(x2)) == f2_24(x2)
    @test f12_24(f2_12(x2)) == f2_24(x2)
    @test f12_24(f4_12(x4)) == f4_24(x4)
    @test f4_12(f2_4(x2)) == f2_12(x2)
    @test f4_24(f2_4(x2)) == f2_24(x2)
    @test f8_16(f4_8(x4)) == f4_16(x4)
    @test f8_24(f4_8(x4)) == f4_24(x4)
    @test f4_8(f2_4(x2)) == f2_8(x2)
    @test f4_16(f2_4(x2)) == f2_16(x2)

    #= Produced an error

    f4_12 = embed(k4, k12)
    f8_16 = embed(k8, k16)
    f2_16 = embed(k2, k16)
    f6_12 = embed(k6, k12) 
    f2_6 = embed(k2, k6)
    f4_16 = embed(k4, k16)
    f2_8 = embed(k2, k8)
    f4_8 = embed(k4, k8)
    f2_12 = embed(k2, k12)
    f2_4 = embed(k2, k4)

    @test f8_16(f2_8(x2)) == f2_16(x2)
    @test f8_16(f4_8(x4)) == f4_16(x4)
    @test f6_12(f2_6(x2)) == f2_12(x2)
    @test f4_12(f2_4(x2)) == f2_12(x2)
    @test f4_8(f2_4(x2)) == f2_8(x2)
    @test f4_16(f2_4(x2)) == f2_16(x2)

    =#

    #=
    @time f35_70 = embed(k35, k70)
    @time f10_70 = embed(k10, k70)
    @time f5_10 = embed(k5, k10)
    @time f5_35 = embed(k5, k35)
    @time f10_210 = embed(k10, k210)
    @time f70_210 = embed(k70, k210)
    @time f5_50 = embed(k5, k50)
    @time f10_50 = embed(k10, k50)

    @test f10_70(f5_10(x5^3)) == f35_70(f5_35(x5^3))
    @test f70_210(f10_70(x10)) == f10_210(x10)
    @test f10_50(f5_10(x5^2+1)) == f5_50(x5^2+1)
    =#

    println("PASS")
end
