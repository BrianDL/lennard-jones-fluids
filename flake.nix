{
  description = "A development environment for Python";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs, ... }@inputs:
  let
    system = "x86_64-linux";
    pkgs = nixpkgs.legacyPackages.${system};
    
    # Define packages as an attribute set
    packageSet = with pkgs; {
      inherit entr;
      testloop = writeShellScriptBin "loop" ./scripts/start_test_loop;
    };

  in rec {
    defaultShell = pkgs.mkShell {
      buildInputs = with packageSet; [
        entr
        testloop
      ];

      shellHook = ''
      '';

    };
    
    ### defining available shells
    devShells.${system} = {
      default = defaultShell;

      allPackages = pkgs.mkShell {
        # Convert the attribute set to a list for buildInputs
        buildInputs = builtins.attrValues packageSet;
      };
    };
    
    ### Assign the attribute set directly to packages
    packages.${system}.default = packageSet;

  };
}