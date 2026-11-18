{ pkgs, lib, config, inputs, ... }:

lib.mkMerge [
  {
    # https://devenv.sh/basics/
    env.LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath [
      pkgs.stdenv.cc.cc.lib
      pkgs.zlib
    ];
    env.PIXI_HOME = "${config.devenv.state}/pixi";

    # https://devenv.sh/packages/
    packages = with pkgs; [
      gh
      pixi
      python312
      python313
      python314
      python315
      python3Packages.tox
    ];

    # https://devenv.sh/languages/
    languages.python = {
      enable = true;
      venv.enable = true;
      manylinux.enable = true;
    };

    # https://devenv.sh/basics/
    enterShell = ''
      if [ -f "pixi.toml" ] || [ -f "pyproject.toml" ]; then
        pixi install
        eval "$(pixi shell-hook)"
      fi
    '';

    # https://devenv.sh/tasks/

    # https://devenv.sh/tests/
    enterTest = ''
      echo "Running tests"
      git --version | grep --color=auto "${pkgs.git.version}"
    '';

    # https://devenv.sh/git-hooks/
    # git-hooks.hooks.shellcheck.enable = true;

    # See full reference at https://devenv.sh/reference/options/
  }
  (lib.mkIf pkgs.stdenv.hostPlatform.isLinux {
    enterShell = ''
      # add libnss_sss.so.2 to LD_PRELOAD if it exists, to avoid issues with SSSD and NSS
      LIBNSS_SSS_PATH=$(find /lib /usr/lib /lib64 /usr/lib64 -name "libnss_sss.so.2" 2>/dev/null | head -n 1)
      if [ -n "$LIBNSS_SSS_PATH" ]; then
        export LD_PRELOAD="$LIBNSS_SSS_PATH"
      fi
    '';
  })
]
