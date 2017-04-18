<?xml version='1.0'?>
<Project Type="Project" LVVersion="8208002">
	<Property Name="NI.Project.Description" Type="Str">LAC demo program that shows how to control and get responses.

Written by: Complete Automated Solutions
www.CompleteAutomatedSolutions.com
</Property>
	<Item Name="My Computer" Type="My Computer">
		<Property Name="server.app.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.control.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.tcp.enabled" Type="Bool">false</Property>
		<Property Name="server.tcp.port" Type="Int">0</Property>
		<Property Name="server.tcp.serviceName" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.tcp.serviceName.default" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.vi.callsEnabled" Type="Bool">true</Property>
		<Property Name="server.vi.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="specify.custom.address" Type="Bool">false</Property>
		<Item Name="docs" Type="Folder">
			<Item Name="Firgelli_LAC_LabVIEW.rtf" Type="Document" URL="docs/Firgelli_LAC_LabVIEW.rtf"/>
			<Item Name="LAC_Advanced_Configuration.pdf" Type="Document" URL="docs/LAC_Advanced_Configuration.pdf"/>
			<Item Name="LAC_datasheet.pdf" Type="Document" URL="docs/LAC_datasheet.pdf"/>
			<Item Name="mpusbapi.h" Type="Document" URL="docs/mpusbapi.h"/>
		</Item>
		<Item Name="USB_VIs" Type="Folder">
			<Item Name="mpusbapi.dll" Type="Document" URL="USB_VIs/mpusbapi.dll"/>
			<Item Name="MPUSBClose.vi" Type="VI" URL="USB_VIs/MPUSBClose.vi"/>
			<Item Name="MPUSBGetDeviceCount.vi" Type="VI" URL="USB_VIs/MPUSBGetDeviceCount.vi"/>
			<Item Name="MPUSBGetDeviceDescriptor.vi" Type="VI" URL="USB_VIs/MPUSBGetDeviceDescriptor.vi"/>
			<Item Name="MPUSBGetDLLVersion.vi" Type="VI" URL="USB_VIs/MPUSBGetDLLVersion.vi"/>
			<Item Name="MPUSBOpen.vi" Type="VI" URL="USB_VIs/MPUSBOpen.vi"/>
			<Item Name="MPUSBRead.vi" Type="VI" URL="USB_VIs/MPUSBRead.vi"/>
			<Item Name="MPUSBWrite.vi" Type="VI" URL="USB_VIs/MPUSBWrite.vi"/>
		</Item>
		<Item Name="supportVIs" Type="Folder">
			<Item Name="actuatorSelect.vi" Type="VI" URL="supportVIs/actuatorSelect.vi"/>
			<Item Name="convertPositionResponse.vi" Type="VI" URL="supportVIs/convertPositionResponse.vi"/>
			<Item Name="convertReadResponse.vi" Type="VI" URL="supportVIs/convertReadResponse.vi"/>
			<Item Name="readPosistionFirgelli.vi" Type="VI" URL="supportVIs/readPosistionFirgelli.vi"/>
			<Item Name="writeAccuracyFirgelli.vi" Type="VI" URL="supportVIs/writeAccuracyFirgelli.vi"/>
			<Item Name="writeLimitsFirgelli.vi" Type="VI" URL="supportVIs/writeLimitsFirgelli.vi"/>
			<Item Name="writePositionFirgelli.vi" Type="VI" URL="supportVIs/writePositionFirgelli.vi"/>
			<Item Name="writeVelocityFirgelli.vi" Type="VI" URL="supportVIs/writeVelocityFirgelli.vi"/>
		</Item>
		<Item Name="Examples" Type="Folder">
			<Item Name="Firgelli_Simple_Configuration.vi" Type="VI" URL="Examples/Firgelli_Simple_Configuration.vi"/>
			<Item Name="Firgelli_Simple_Limits.vi" Type="VI" URL="Examples/Firgelli_Simple_Limits.vi"/>
			<Item Name="Firgelli_Simple_Move.vi" Type="VI" URL="Examples/Firgelli_Simple_Move.vi"/>
		</Item>
		<Item Name="Firgelli_LAC_LabVIEW.vi" Type="VI" URL="Firgelli_LAC_LabVIEW.vi"/>
		<Item Name="app.ico" Type="Document" URL="app.ico"/>
		<Item Name="Dependencies" Type="Dependencies">
			<Item Name="vi.lib" Type="Folder">
				<Item Name="General Error Handler.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/General Error Handler.vi"/>
				<Item Name="General Error Handler CORE.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/General Error Handler CORE.vi"/>
				<Item Name="Check Special Tags.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Check Special Tags.vi"/>
				<Item Name="TagReturnType.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/TagReturnType.ctl"/>
				<Item Name="Set String Value.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Set String Value.vi"/>
				<Item Name="GetRTHostConnectedProp.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/GetRTHostConnectedProp.vi"/>
				<Item Name="Error Code Database.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Error Code Database.vi"/>
				<Item Name="Trim Whitespace.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Trim Whitespace.vi"/>
				<Item Name="whitespace.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/whitespace.ctl"/>
				<Item Name="Format Message String.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Format Message String.vi"/>
				<Item Name="Find Tag.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Find Tag.vi"/>
				<Item Name="Search and Replace Pattern.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Search and Replace Pattern.vi"/>
				<Item Name="Set Bold Text.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Set Bold Text.vi"/>
				<Item Name="Details Display Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Details Display Dialog.vi"/>
				<Item Name="Clear Errors.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Clear Errors.vi"/>
				<Item Name="DialogTypeEnum.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/DialogTypeEnum.ctl"/>
				<Item Name="ErrWarn.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/ErrWarn.ctl"/>
				<Item Name="eventvkey.ctl" Type="VI" URL="/&lt;vilib&gt;/event_ctls.llb/eventvkey.ctl"/>
				<Item Name="Not Found Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Not Found Dialog.vi"/>
				<Item Name="Three Button Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Three Button Dialog.vi"/>
				<Item Name="Three Button Dialog CORE.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Three Button Dialog CORE.vi"/>
				<Item Name="Longest Line Length in Pixels.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Longest Line Length in Pixels.vi"/>
				<Item Name="Convert property node font to graphics font.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Convert property node font to graphics font.vi"/>
				<Item Name="Get Text Rect.vi" Type="VI" URL="/&lt;vilib&gt;/picture/picture.llb/Get Text Rect.vi"/>
				<Item Name="Get String Text Bounds.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Get String Text Bounds.vi"/>
				<Item Name="LVBoundsTypeDef.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/miscctls.llb/LVBoundsTypeDef.ctl"/>
				<Item Name="BuildHelpPath.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/BuildHelpPath.vi"/>
				<Item Name="GetHelpDir.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/GetHelpDir.vi"/>
				<Item Name="DialogType.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/DialogType.ctl"/>
			</Item>
			<Item Name="controlValues.ctl" Type="VI" URL="typeDefs/controlValues.ctl"/>
		</Item>
		<Item Name="Build Specifications" Type="Build">
			<Item Name="Firgelli_LAC_LabVIEW" Type="EXE">
				<Property Name="App_copyErrors" Type="Bool">true</Property>
				<Property Name="App_INI_aliasGUID" Type="Str">{BF701982-0E2C-4455-AD07-7D72776D9399}</Property>
				<Property Name="App_INI_GUID" Type="Str">{E2EA0623-80B8-4FCE-80E2-59DC343D1B46}</Property>
				<Property Name="Bld_buildSpecName" Type="Str">Firgelli_LAC_LabVIEW</Property>
				<Property Name="Bld_excludeLibraryItems" Type="Bool">true</Property>
				<Property Name="Bld_excludePolymorphicVIs" Type="Bool">true</Property>
				<Property Name="Bld_localDestDir" Type="Path">../builds/Firgelli_LAC_LabVIEW/Firgelli_LAC_LabVIEW</Property>
				<Property Name="Bld_localDestDirType" Type="Str">relativeToCommon</Property>
				<Property Name="Bld_modifyLibraryFile" Type="Bool">true</Property>
				<Property Name="Destination[0].destName" Type="Str">Firgelli_LAC_LabVIEW.exe</Property>
				<Property Name="Destination[0].path" Type="Path">../builds/NI_AB_PROJECTNAME/Firgelli_LAC_LabVIEW/Firgelli_LAC_LabVIEW.exe</Property>
				<Property Name="Destination[0].preserveHierarchy" Type="Bool">true</Property>
				<Property Name="Destination[0].type" Type="Str">App</Property>
				<Property Name="Destination[1].destName" Type="Str">Support Directory</Property>
				<Property Name="Destination[1].path" Type="Path">../builds/NI_AB_PROJECTNAME/Firgelli_LAC_LabVIEW/data</Property>
				<Property Name="DestinationCount" Type="Int">2</Property>
				<Property Name="Exe_iconItemID" Type="Ref">/My Computer/app.ico</Property>
				<Property Name="Source[0].itemID" Type="Str">{2236AC68-2938-4216-BC2F-B3879618B75F}</Property>
				<Property Name="Source[0].type" Type="Str">Container</Property>
				<Property Name="Source[1].destinationIndex" Type="Int">0</Property>
				<Property Name="Source[1].itemID" Type="Ref">/My Computer/Firgelli_LAC_LabVIEW.vi</Property>
				<Property Name="Source[1].sourceInclusion" Type="Str">TopLevel</Property>
				<Property Name="Source[1].type" Type="Str">VI</Property>
				<Property Name="Source[2].destinationIndex" Type="Int">0</Property>
				<Property Name="Source[2].itemID" Type="Ref">/My Computer/USB_VIs/mpusbapi.dll</Property>
				<Property Name="Source[2].sourceInclusion" Type="Str">Include</Property>
				<Property Name="SourceCount" Type="Int">3</Property>
				<Property Name="TgtF_companyName" Type="Str">CAS</Property>
				<Property Name="TgtF_fileDescription" Type="Str">Firgelli_LAC_LabVIEW</Property>
				<Property Name="TgtF_fileVersion.major" Type="Int">1</Property>
				<Property Name="TgtF_internalName" Type="Str">Firgelli_LAC_LabVIEW</Property>
				<Property Name="TgtF_legalCopyright" Type="Str">Copyright © 2011 CAS</Property>
				<Property Name="TgtF_productName" Type="Str">Firgelli_LAC_LabVIEW</Property>
				<Property Name="TgtF_targetfileGUID" Type="Str">{F55076E0-1832-4252-88E1-96FA05A03692}</Property>
				<Property Name="TgtF_targetfileName" Type="Str">Firgelli_LAC_LabVIEW.exe</Property>
			</Item>
			<Item Name="Firgelli_LAC_LabVIEW Installer" Type="Installer">
				<Property Name="Destination[0].name" Type="Str">Firgelli_LAC_LabVIEW</Property>
				<Property Name="Destination[0].parent" Type="Str">{3912416A-D2E5-411B-AFEE-B63654D690C0}</Property>
				<Property Name="Destination[0].tag" Type="Str">{1053A83D-16C5-4919-9B27-B0203BEA59D3}</Property>
				<Property Name="Destination[0].type" Type="Str">userFolder</Property>
				<Property Name="DestinationCount" Type="Int">1</Property>
				<Property Name="INST_author" Type="Str">CAS</Property>
				<Property Name="INST_buildLocation" Type="Path">../builds/Firgelli_LAC_LabVIEW/Firgelli_LAC_LabVIEW Installer</Property>
				<Property Name="INST_buildLocation.type" Type="Str">relativeToCommon</Property>
				<Property Name="INST_buildSpecName" Type="Str">Firgelli_LAC_LabVIEW Installer</Property>
				<Property Name="INST_defaultDir" Type="Str">{1053A83D-16C5-4919-9B27-B0203BEA59D3}</Property>
				<Property Name="INST_productName" Type="Str">Firgelli_LAC_LabVIEW</Property>
				<Property Name="INST_productVersion" Type="Str">1.0.0</Property>
				<Property Name="InstSpecBitness" Type="Str">32-bit</Property>
				<Property Name="InstSpecVersion" Type="Str">10018002</Property>
				<Property Name="MSI_arpCompany" Type="Str">Complete Automated Solutions</Property>
				<Property Name="MSI_arpContact" Type="Str">Matt Fitzsimons</Property>
				<Property Name="MSI_arpPhone" Type="Str">847-828-3389</Property>
				<Property Name="MSI_arpURL" Type="Str">http://www.completeAutomatedSolutions.com/</Property>
				<Property Name="MSI_distID" Type="Str">{A8BAEB53-3102-4D4B-A4C5-0F7ACA4A10C2}</Property>
				<Property Name="MSI_osCheck" Type="Int">0</Property>
				<Property Name="MSI_upgradeCode" Type="Str">{E3F93C50-7081-40FE-86FA-AE812C1C547E}</Property>
				<Property Name="MSI_windowMessage" Type="Str">If you run the executable demo program on a system that LabVIEW 2010 has not been installed you will need to install:

LabVIEW Run-Time Engine 2010 - (32-bit Standard RTE) - Windows 7/7 64 bit/Server 2003 R2 (32-bit)/XP/Vista/Vista x64/Server 2008 R2 (64-bit)

http://joule.ni.com/nidu/cds/view/p/id/2087/lang/en

Example Program Written by:

Complete Automated Solutions
http://www.completeautomatedsolutions.com</Property>
				<Property Name="MSI_windowTitle" Type="Str">Firgelli LabVIEW Sample Program</Property>
				<Property Name="RegDest[0].dirName" Type="Str">Software</Property>
				<Property Name="RegDest[0].dirTag" Type="Str">{DDFAFC8B-E728-4AC8-96DE-B920EBB97A86}</Property>
				<Property Name="RegDest[0].parentTag" Type="Str">2</Property>
				<Property Name="RegDestCount" Type="Int">1</Property>
				<Property Name="Source[0].dest" Type="Str">{1053A83D-16C5-4919-9B27-B0203BEA59D3}</Property>
				<Property Name="Source[0].File[0].dest" Type="Str">{1053A83D-16C5-4919-9B27-B0203BEA59D3}</Property>
				<Property Name="Source[0].File[0].name" Type="Str">Firgelli_LAC_LabVIEW.exe</Property>
				<Property Name="Source[0].File[0].Shortcut[0].destIndex" Type="Int">0</Property>
				<Property Name="Source[0].File[0].Shortcut[0].name" Type="Str">Firgelli_LAC_LabVIEW</Property>
				<Property Name="Source[0].File[0].Shortcut[0].subDir" Type="Str"></Property>
				<Property Name="Source[0].File[0].ShortcutCount" Type="Int">1</Property>
				<Property Name="Source[0].File[0].tag" Type="Str">{F55076E0-1832-4252-88E1-96FA05A03692}</Property>
				<Property Name="Source[0].FileCount" Type="Int">1</Property>
				<Property Name="Source[0].name" Type="Str">Firgelli_LAC_LabVIEW</Property>
				<Property Name="Source[0].tag" Type="Ref">/My Computer/Build Specifications/Firgelli_LAC_LabVIEW</Property>
				<Property Name="Source[0].type" Type="Str">EXE</Property>
				<Property Name="SourceCount" Type="Int">1</Property>
			</Item>
		</Item>
	</Item>
</Project>
