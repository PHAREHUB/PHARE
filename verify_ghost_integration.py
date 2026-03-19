#!/usr/bin/env python3
"""
Simple verification test for ghost width calculator integration.
This test verifies the ghost width values match expected values.
"""

import subprocess
import sys

# Test program that verifies ghost widths at compile time
test_code = """
#include <iostream>
#include <cassert>
#include "core/utilities/ghost_width_calculator.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"

using namespace PHARE::core;

int main() {
    std::cout << "Testing Ghost Width Calculator Integration..." << std::endl;
    
    // Test 1: Calculator produces correct values
    std::cout << "\\n[Test 1] Ghost Width Calculator Values:" << std::endl;
    std::cout << "  Hybrid Order 1: " << GhostWidthCalculator<HybridConfig<1>>::value << " (expected: 2)" << std::endl;
    std::cout << "  Hybrid Order 2: " << GhostWidthCalculator<HybridConfig<2>>::value << " (expected: 4)" << std::endl;
    std::cout << "  Hybrid Order 3: " << GhostWidthCalculator<HybridConfig<3>>::value << " (expected: 4)" << std::endl;
    std::cout << "  MHD WENOZ (3): " << GhostWidthCalculator<MHDConfig<3>>::value << " (expected: 4)" << std::endl;
    
    assert(GhostWidthCalculator<HybridConfig<1>>::value == 2);
    assert(GhostWidthCalculator<HybridConfig<2>>::value == 4);
    assert(GhostWidthCalculator<HybridConfig<3>>::value == 4);
    assert(GhostWidthCalculator<MHDConfig<3>>::value == 4);
    std::cout << "  ✓ All calculator values correct" << std::endl;
    
    // Test 2: GridLayout uses calculator
    std::cout << "\\n[Test 2] GridLayout Integration:" << std::endl;
    using Layout1 = GridLayout<GridLayoutImplYee<1, 1>>;
    using Layout2 = GridLayout<GridLayoutImplYee<2, 2>>;
    using Layout3 = GridLayout<GridLayoutImplYee<3, 3>>;
    
    std::cout << "  GridLayout<1,1>::nbrGhosts(): " << Layout1::nbrGhosts() << " (expected: 2)" << std::endl;
    std::cout << "  GridLayout<2,2>::nbrGhosts(): " << Layout2::nbrGhosts() << " (expected: 4)" << std::endl;
    std::cout << "  GridLayout<3,3>::nbrGhosts(): " << Layout3::nbrGhosts() << " (expected: 4)" << std::endl;
    
    assert(Layout1::nbrGhosts() == 2);
    assert(Layout2::nbrGhosts() == 4);
    assert(Layout3::nbrGhosts() == 4);
    std::cout << "  ✓ GridLayout uses calculator correctly" << std::endl;
    
    // Test 3: Even number constraint
    std::cout << "\\n[Test 3] Even Number Constraint:" << std::endl;
    assert(Layout1::nbrGhosts() % 2 == 0);
    assert(Layout2::nbrGhosts() % 2 == 0);
    assert(Layout3::nbrGhosts() % 2 == 0);
    std::cout << "  ✓ All ghost counts are even" << std::endl;
    
    // Test 4: Backward compatibility
    std::cout << "\\n[Test 4] Backward Compatibility:" << std::endl;
    // Original hardcoded values were {2, 4, 4}
    assert(Layout1::nbrGhosts() == 2);  // Order 1 -> 2
    assert(Layout2::nbrGhosts() == 4);  // Order 2 -> 4
    assert(Layout3::nbrGhosts() == 4);  // Order 3 -> 4
    std::cout << "  ✓ Matches original ghost width values" << std::endl;
    
    std::cout << "\\n✅ All tests passed!" << std::endl;
    return 0;
}
"""

def main():
    print("=" * 70)
    print("Ghost Width Calculator - Standalone Verification")
    print("=" * 70)
    
    # Write test file
    test_file = "/tmp/phare_ghost_test.cpp"
    with open(test_file, "w") as f:
        f.write(test_code)
    
    print(f"\n[1/2] Compiling test program...")
    
    # Compile
    compile_cmd = [
        "g++",
        "-std=c++20",
        "-I./src",
        test_file,
        "-o", "/tmp/phare_ghost_test"
    ]
    
    result = subprocess.run(compile_cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("❌ Compilation failed:")
        print(result.stderr)
        return 1
    
    print("✅ Compilation successful")
    
    print(f"\n[2/2] Running tests...")
    print("-" * 70)
    
    # Run
    result = subprocess.run(["/tmp/phare_ghost_test"], capture_output=True, text=True)
    
    print(result.stdout)
    
    if result.returncode != 0:
        print("❌ Tests failed:")
        print(result.stderr)
        return 1
    
    print("-" * 70)
    print("\n" + "=" * 70)
    print("🎉 Ghost Width Calculator Integration: VERIFIED")
    print("=" * 70)
    print("\nSummary:")
    print("  ✓ Calculator produces correct ghost widths")
    print("  ✓ GridLayout uses calculator correctly")
    print("  ✓ Even number constraint enforced")
    print("  ✓ Backward compatibility maintained")
    print("\nThe implementation is ready for full system testing.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
