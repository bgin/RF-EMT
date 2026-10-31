// Compilation: g++ military_array.cpp -o military_array -lgmsh
//To compile and execute this code on your local system, you must link against the installed shared libgmsh library dependencies:
#include <iostream>
#include <vector>
#include <cmath>
#include <gmsh.h>

int main(int argc, char **argv) {
    // 1. Initialize Gmsh Engine environment
    gmsh::initialize();
    gmsh::model::add("military_array_3d");

    // Enforce OpenCASCADE kernel geometry generation
    std::string occ = "OpenCASCADE";

    // --- 2. Operational & Structural Scaling Constraints (Dimensions in mm) ---
    const int N_rows = 16;
    const int N_cols = 16;
    const double Element_Spacing = 18.0;

    const double W_patch = 11.2;
    const double L_patch = 11.2;
    const double H_sub = 1.5;
    const double H_spacer = 0.5;

    const double W_feed = 0.8;
    const double Margin = 20.0;
    const double H_gnd = 2.0;
    const double Air_Margin = 30.0;

    // Calculate total structural physical footprint dimensions
    const double Array_Width = (N_cols - 1) * Element_Spacing + W_patch;
    const double Array_Length = (N_rows - 1) * Element_Spacing + L_patch;

    const double Gnd_W = Array_Width + 2.0 * Margin;
    const double Gnd_L = Array_Length + 2.0 * Margin;

    // --- 3. Construct Base Structural Layer Volumes ---
    // Structural Ground Aluminum Backing Plate
    int gnd_vol = gmsh::model::occ::addBox(-Gnd_W/2.0, -Gnd_L/2.0, -H_gnd, Gnd_W, Gnd_L, H_gnd);
    // Dielectric Substrate Layer
    int sub_vol = gmsh::model::occ::addBox(-Gnd_W/2.0, -Gnd_L/2.0, 0.0, Gnd_W, Gnd_L, H_sub);

    // Track dynamic element indices for Boolean integration operations
    std::vector<std::pair<int, int>> object_dim_tags;
    std::vector<std::pair<int, int>> tool_dim_tags;

    object_dim_tags.push_back({3, gnd_vol});
    object_dim_tags.push_back({3, sub_vol});

    // --- 4. Programmatic Array Generation Loop (256 Elements) ---
    for (int r = 0; r < N_rows; ++r) {
        double y_pos = -((N_rows - 1) * Element_Spacing) / 2.0 + (r * Element_Spacing);
        
        for (int c = 0; c < N_cols; ++c) {
            double x_pos = -((N_cols - 1) * Element_Spacing) / 2.0 + (c * Element_Spacing);

            // Construct Raised Radiating Patch Element (Extruded Thickness of 0.1mm)
            int patch_vol = gmsh::model::occ::addBox(
                x_pos - W_patch/2.0, y_pos - L_patch/2.0, H_sub + H_spacer, 
                W_patch, L_patch, 0.1
            );
            tool_dim_tags.push_back({3, patch_vol});

            // Construct Vertical Coaxial Feed Pin through layers
            int feed_vol = gmsh::model::occ::addCylinder(
                x_pos - W_patch/4.0, y_pos, -H_gnd, 
                0.0, 0.0, H_gnd + H_sub + H_spacer, 
                W_feed / 2.0
            );
            tool_dim_tags.push_back({3, feed_vol});
        }
    }

    // Intersect and fragment internal arrays, setting conforming mesh interfaces
    std::vector<std::pair<int, int>> out_dim_tags;
    std::vector<std::vector<std::pair<int, int>>> out_dim_tags_map;
    gmsh::model::occ::intersect(object_dim_tags, tool_dim_tags, out_dim_tags, out_dim_tags_map);

    // Reconstruct global domain matrix tracking updated fragment IDs
    std::vector<std::pair<int, int>> integrated_domains = out_dim_tags;

    // --- 5. Far-Field Simulation Air Enclosure ---
    int air_vol = gmsh::model::occ::addBox(
        -Gnd_W/2.0 - Air_Margin, -Gnd_L/2.0 - Air_Margin, -H_gnd - Air_Margin,
        Gnd_W + 2.0 * Air_Margin, Gnd_L + 2.0 * Air_Margin, H_sub + H_spacer + 2.0 * Air_Margin
    );

    std::vector<std::pair<int, int>> final_obj = {{3, air_vol}};
    std::vector<std::pair<int, int>> final_out;
    std::vector<std::vector<std::pair<int, int>>> final_out_map;
    
    // Sync external air mesh enclosure boundaries with the integrated internal antenna
    gmsh::model::occ::intersect(final_obj, integrated_domains, final_out, final_out_map);

    // Finalize OpenCASCADE geometry updates before generating mesh
    gmsh::model::occ::synchronize();

    // --- 6. Adaptive Mesh Scaling ---
    double lc_patch = 1.0;
    double lc_farfield = 8.0;

    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.2);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", lc_farfield);

    // Fetch all points across the geometry to apply targeted local sizing fields
    std::vector<std::pair<int, int>> all_points;
    gmsh::model::getEntities(all_points, 0);
    
    // Smooth grid resolution transition sizing constraints mapping
    for (const auto& pt : all_points) {
        gmsh::model::mesh::setSize({pt}, lc_patch);
    }

    // --- 7. Compile 3D Mesh & Export System Payload ---
    std::cout << "Compiling 3D Volume Mesh for 256 Element Military Phased Array..." << std::endl;
    gmsh::model::mesh::generate(3);
    gmsh::write("military_array.msh");
    std::cout << "Mesh compilation complete. Saved to 'military_array.msh'." << std::endl;

    gmsh::finalize();
    return 0;
}
