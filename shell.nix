{
  pkgs ? import <nixpkgs> {},
}:

let
  # Use builtins.getFlake to obtain the local flake's outputs
  flake = builtins.getFlake (toString ./.);

  # Access the packageSet defined in the flake
  devShell = flake.outputs.devShells.x86_64-linux.default;
  # packageSet = flake.outputs.packages.${flake.inputs.nixpkgs.system}.packageSet;

in devShell