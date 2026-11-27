
#include <cstdint>
#include <cstdlib>

#include <core/Environment.h>
#include <core/debug/TestSubsystem.h>
#include <core/DataTypes.h>

#include <field/Field.h>
#include <field/GhostLayerField.h>
#include <field/Layout.h>

#include "walberla_helper/field/Field.h"

namespace field_test {

    using namespace walberla;

    const field::Layout layout{ field::Layout::fzyx };

    const cell_idx_t xs = 2;
    const cell_idx_t ys = 4;
    const cell_idx_t zs = 8;
    const cell_idx_t fs = 3;

    void compareAgainstFieldTest() {

        Field<int, fs> standardField(xs, ys, zs, layout);

        // initialise field with random data
        for( auto it = standardField.begin(); it != standardField.end(); ++it ) {
            *it = std::rand();
        }

        // extract field info to custom field
        turbine_core::field::Field<int> customField(&standardField);

        // check for equality of data
        for( cell_idx_t f = 0; f < fs; ++f ) {
            for( cell_idx_t z = 0; z < zs; ++z ) {
                for( cell_idx_t y = 0; y < ys; ++y ) {
                    for( cell_idx_t x = 0; x < xs; ++x ) {
                        WALBERLA_CHECK_EQUAL( standardField.get(x,y,z,f), customField.get(x,y,z,f) );
                    }
                }
            }
        }
    }

    void compareAgainstGhostFieldTest( const uint_t nGhostLayers ) {

        GhostLayerField<int, fs> standardField(xs, ys, zs, nGhostLayers, 0, layout);

        // initialise field with random data
        for( auto it = standardField.begin(); it != standardField.end(); ++it ) {
            *it = std::rand();
        }

        // extract field info to custom field
        turbine_core::field::Field<int> customField(&standardField);

        // check for equality of data
        for( cell_idx_t f = 0; f < fs; ++f ) {
            for( cell_idx_t z = 0; z < zs; ++z ) {
                for( cell_idx_t y = 0; y < ys; ++y ) {
                    for( cell_idx_t x = 0; x < xs; ++x ) {
                        WALBERLA_CHECK_EQUAL( standardField.get(x,y,z,f), customField.get(x,y,z,f) );
                    }
                }
            }
        }

    }

    template< typename Type_T >
    struct IndexOnlyField {

        IndexOnlyField( uint_t xAlloc, uint_t yAlloc, uint_t zAlloc, uint_t f, cell_idx_t xo, cell_idx_t yo, cell_idx_t zo )
            : xA_(xAlloc), yA_(yAlloc), zA_(zAlloc), f_(f), xo_(xo), yo_(yo), zo_(zo) {}

        uint_t xSize() const { return xA_ - uint_t(2*xo_); }
        uint_t ySize() const { return yA_ - uint_t(2*yo_); }
        uint_t zSize() const { return zA_ - uint_t(2*zo_); }
        uint_t fSize() const { return f_; }
        uint_t xAllocSize() const { return xA_; }
        uint_t yAllocSize() const { return yA_; }
        uint_t zAllocSize() const { return zA_; }
        uint_t fAllocSize() const { return f_; }
        int64_t xStride() const { return 1; }
        int64_t yStride() const { return int64_t(xA_); }
        int64_t zStride() const { return int64_t(xA_) * int64_t(yA_); }
        int64_t fStride() const { return int64_t(xA_) * int64_t(yA_) * int64_t(zA_); }
        cell_idx_t xOff() const { return xo_; }
        cell_idx_t yOff() const { return yo_; }
        cell_idx_t zOff() const { return zo_; }
        field::Layout layout() const { return field::Layout::fzyx; }
        Type_T * dataAt( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) const {
            return base() + ( int64_t(f) * fStride() + int64_t(z) * zStride() + int64_t(y) * yStride() + int64_t(x) * xStride() );
        }

        static Type_T * base() { return reinterpret_cast<Type_T*>( std::uintptr_t(1) << 40 ); }

    private:
        uint_t xA_, yA_, zA_, f_;
        cell_idx_t xo_, yo_, zo_;
    };

    void indexArithmeticTest( uint_t xAlloc, uint_t yAlloc, uint_t zAlloc, uint_t fSize,
                              cell_idx_t xo, cell_idx_t yo, cell_idx_t zo ) {

        IndexOnlyField<float> mock( xAlloc, yAlloc, zAlloc, fSize, xo, yo, zo );
        const turbine_core::field::Field<float> customField(&mock);

        const std::int64_t xS = 1, yS = std::int64_t(xAlloc), zS = std::int64_t(xAlloc) * std::int64_t(yAlloc), fS = zS * std::int64_t(zAlloc);

        const cell_idx_t xs = cell_idx_t(xAlloc) - 2*xo, ys = cell_idx_t(yAlloc) - 2*yo, zs = cell_idx_t(zAlloc) - 2*zo;

        const cell_idx_t xList[3] = { -xo, xs / 2, xs - 1 + xo };
        const cell_idx_t yList[3] = { -yo, ys / 2, ys - 1 + yo };
        const cell_idx_t zList[3] = { -zo, zs / 2, zs - 1 + zo };
        const cell_idx_t fList[3] = { 0, cell_idx_t(fSize / 2), cell_idx_t(fSize) - 1 };

        for( cell_idx_t f : fList ) for( cell_idx_t z : zList ) for( cell_idx_t y : yList ) for( cell_idx_t x : xList ) {
            const std::int64_t expected = std::int64_t(f) * fS + std::int64_t(x) * xS + std::int64_t(y) * yS + std::int64_t(z) * zS;
            const float * address = &customField.get(x,y,z,f);
            const std::int64_t offset = std::int64_t( address - IndexOnlyField<float>::base() );
            WALBERLA_CHECK_EQUAL( offset, expected );
        }
    }

    void largeIndexTest() {
        indexArithmeticTest( 442, 442, 442, 27, 1, 1, 1 );
        indexArithmeticTest( 1302, 1302, 1302, 1, 1, 1, 1 );
        indexArithmeticTest( 1304, 1302, 1300, 3, 2, 1, 0 );
    }

    void constGhostLayersTest() {
        Field<int, fs> standardField(xs, ys, zs, layout);
        const turbine_core::field::Field<int> customField(&standardField);
        WALBERLA_CHECK_EQUAL( customField.nrOfGhostLayers(), cell_idx_t(0) );
    }

    //TODO setup test case for GPUField

    int main( int argc, char** argv ) {
        Environment walberlaEnv( argc, argv );

        debug::enterTestMode();

        compareAgainstFieldTest();

        compareAgainstGhostFieldTest(0);
        compareAgainstGhostFieldTest(1);
        compareAgainstGhostFieldTest(2);
        compareAgainstGhostFieldTest(4);

        largeIndexTest();
        constGhostLayersTest();

        return EXIT_SUCCESS;
    }

}

int main( int argc, char** argv ) {
    return field_test::main(argc, argv);
}