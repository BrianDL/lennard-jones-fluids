{
  description = "A development environment for Python";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, ... }@inputs:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        
        # Define packages as an attribute set
        packageSet = with pkgs; {
          inherit entr;
        };

      in rec {
        defaultShell = pkgs.mkShell {
          buildInputs = with packageSet; [
            entr
          ];

        };
        
        ### defining available shells
        devShells = {
          default = defaultShell;

          allPackages = pkgs.mkShell {
            # Convert the attribute set to a list for buildInputs
            buildInputs = builtins.attrValues packageSet;
          };
        };
        
        ### Assign the attribute set directly to packages
        packages = packageSet;

        # defaultPackage = self.packages.${system}.deployRocNightly;
      }
    );
}